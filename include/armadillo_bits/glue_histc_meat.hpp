// Copyright (C) 2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup glue_histc
//! @{


template<typename eT>
inline
void
glue_histc::apply_noalias(Mat<uword>& C, const Mat<eT>& A, const Mat<eT>& B, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ((B.is_vec() == false) && (B.is_empty() == false)), "histc(): parameter 'edges' is not a vector" );
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_elem = B.n_elem;
  
  if( B_n_elem == uword(0) )  { C.reset(); return; }
  
  const eT*   B_mem       = B.memptr();
  const uword B_n_elem_m1 = B_n_elem - 1;
  
  if(dim == uword(0))
    {
    C.zeros(B_n_elem, A_n_cols);
    
    for(uword col=0; col < A_n_cols; ++col)
      {
      const eT*    A_coldata = A.colptr(col);
            uword* C_coldata = C.colptr(col);
      
      for(uword row=0; row < A_n_rows; ++row)
        {
        const eT x = A_coldata[row];
        
        for(uword i=0; i < B_n_elem_m1; ++i)
          {
               if( (B_mem[i]           <= x) && (x < B_mem[i+1]) )  { C_coldata[i]++;           break; }
          else if(  B_mem[B_n_elem_m1] == x                      )  { C_coldata[B_n_elem_m1]++; break; }    // for compatibility with Matlab
          }
        }
      }
    }
  else
  if(dim == uword(1))
    {
    C.zeros(A_n_rows, B_n_elem);
    
    if(A.n_rows == 1)
      {
      const uword  A_n_elem = A.n_elem;
      const eT*    A_mem    = A.memptr();
            uword* C_mem    = C.memptr();
      
      for(uword j=0; j < A_n_elem; ++j)
        {
        const eT x = A_mem[j];
        
        for(uword i=0; i < B_n_elem_m1; ++i)
          {
               if( (B_mem[i]           <= x) && (x < B_mem[i+1]) )  { C_mem[i]++;           break; }
          else if(  B_mem[B_n_elem_m1] == x                      )  { C_mem[B_n_elem_m1]++; break; }    // for compatibility with Matlab
          }
        }
      }
    else
      {
      for(uword row=0; row < A_n_rows; ++row)
      for(uword col=0; col < A_n_cols; ++col)
        {
        const eT x = A.at(row,col);
        
        for(uword i=0; i < B_n_elem_m1; ++i)
          {
               if( (B_mem[i]            <= x) && (x < B_mem[i+1]) )  { C.at(row,i)++;           break; }
          else if(  B_mem[B_n_elem_m1]  == x                      )  { C.at(row,B_n_elem_m1)++; break; }   // for compatibility with Matlab
          }
        }
      }
    }
  }



template<typename T1, typename T2>
inline
void
glue_histc::apply(Mat<uword>& C, const mtGlue<uword,T1,T2,glue_histc>& expr)
  {
  arma_extra_debug_sigprint();
  
  const uword dim = expr.aux_uword;
  
  arma_debug_check( (dim > 1), "histc(): parameter 'dim' must be 0 or 1" );
  
  const quasi_unwrap<T1> UA(expr.A);
  const quasi_unwrap<T2> UB(expr.B);
  
  if(UA.is_alias(C) || UB.is_alias(C))
    {
    Mat<uword> tmp;
    
    glue_histc::apply_noalias(tmp, UA.M, UB.M, dim);
    
    C.steal_mem(tmp);
    }
  else
    {
    glue_histc::apply_noalias(C, UA.M, UB.M, dim);
    }
  }



template<typename T1, typename T2>
inline
void
glue_histc_default::apply(Mat<uword>& C, const mtGlue<uword,T1,T2,glue_histc_default>& expr)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> UA(expr.A);
  const quasi_unwrap<T2> UB(expr.B);
  
  //const uword dim = ( (T1::is_row) || ((UA.M.vec_state == 0) && (UA.M.n_elem <= 1) && (C.vec_state == 2)) ) ? 1 : 0;
  const uword dim = (T1::is_row) ? 1 : 0;
  
  if(UA.is_alias(C) || UB.is_alias(C))
    {
    Mat<uword> tmp;
    
    glue_histc::apply_noalias(tmp, UA.M, UB.M, dim);
    
    C.steal_mem(tmp);
    }
  else
    {
    glue_histc::apply_noalias(C, UA.M, UB.M, dim);
    }
  }


//! @}
