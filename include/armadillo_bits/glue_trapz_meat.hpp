// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup glue_trapz
//! @{



template<typename T1, typename T2>
inline
void
glue_trapz::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_trapz>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword;
  
  const unwrap<T1> UX(in.A);
  const unwrap<T2> UY(in.B);
  
  const Mat<eT>& X = UX.M;
  const Mat<eT>& Y = UY.M;
  
  if( (&out == &X) || (&out == &Y) )
    {
    Mat<eT> tmp;
    
    glue_trapz::apply_noalias(tmp, X, Y, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_trapz::apply_noalias(out, X, Y, dim);
    }
  }



template<typename eT>
inline
void
glue_trapz::apply_noalias(Mat<eT>& out, const Mat<eT>& X, const Mat<eT>& Y, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (dim > 1), "trapz(): argument 'dim' must be 0 or 1" );
  
  arma_debug_check( ((X.is_vec() == false) && (X.is_empty() == false)), "trapz(): argument 'X' must be a vector" );
  
  const uword N = X.n_elem;
  
  if(dim == 0)
    {
    arma_debug_check( (N != Y.n_rows), "trapz(): length of X must equal the number of rows in Y when dim=0" );
    }
  else
  if(dim == 1)
    {
    arma_debug_check( (N != Y.n_cols), "trapz(): length of X must equal the number of columns in Y when dim=1" );
    }
  
  if(N <= 1)
    {
         if(dim == 0)  { out.zeros(1, Y.n_cols); }
    else if(dim == 1)  { out.zeros(Y.n_rows, 1); }
    
    return;
    }
  
  const Col<eT> vec_X( const_cast<eT*>(X.memptr()), X.n_elem, false, true );
  
  const Col<eT> diff_X = diff(vec_X);
  
  if(dim == 0)
    {
    const Row<eT> diff_X_t( const_cast<eT*>(diff_X.memptr()), diff_X.n_elem, false, true );
    
    out = diff_X_t * (0.5 * (Y.rows(0, N-2) + Y.rows(1, N-1)));
    }
  else
  if(dim == 1)
    {
    out = (0.5 * (Y.cols(0, N-2) + Y.cols(1, N-1))) * diff_X;
    }
  }



template<typename T1>
inline
void
op_trapz::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trapz>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  
  const unwrap<T1>   UY(in.m);
  const Mat<eT>& Y = UY.M;
  
  if(&out == &Y)
    {
    Mat<eT> tmp;
    
    op_trapz::apply_noalias(tmp, Y, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_trapz::apply_noalias(out, Y, dim);
    }
  }



template<typename eT>
inline
void
op_trapz::apply_noalias(Mat<eT>& out, const Mat<eT>& Y, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (dim > 1), "trapz(): argument 'dim' must be 0 or 1" );
  
  uword N = 0;
  
       if(dim == 0)  { N = Y.n_rows; }
  else if(dim == 1)  { N = Y.n_cols; }
  
  if(N <= 1)
    {
         if(dim == 0)  { out.zeros(1, Y.n_cols); }
    else if(dim == 1)  { out.zeros(Y.n_rows, 1); }
    
    return;
    }
  
  if(dim == 0)
    {
    out = sum( (0.5 * (Y.rows(0, N-2) + Y.rows(1, N-1))), 0 );
    }
  else
  if(dim == 1)
    {
    out = sum( (0.5 * (Y.cols(0, N-2) + Y.cols(1, N-1))), 1 );
    }
  }



//! @}
