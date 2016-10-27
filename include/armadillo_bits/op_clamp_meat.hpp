// Copyright (C) 2014-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup op_clamp
//! @{



template<typename T1>
inline
void
op_clamp::apply(Mat<typename T1::elem_type>& out, const mtOp<typename T1::elem_type, T1, op_clamp>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(in.m);
  
  if(P.is_alias(out) && (is_Mat<T1>::value == false))
    {
    Mat<eT> tmp;
    
    op_clamp::apply_noalias(tmp, P, in.aux, in.aux_out_eT);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_clamp::apply_noalias(out, P, in.aux, in.aux_out_eT);
    }
  }



template<typename T1>
inline
void
op_clamp::apply_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const typename T1::elem_type min_val, const typename T1::elem_type max_val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
  eT* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    const uword N = P.get_n_elem();
    
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    uword j;
    for(j=1; j<N; j+=2)
      {
      eT val_i = A[j-1];
      eT val_j = A[j  ];
      
      val_i = (val_i < min_val) ? min_val : ((val_i > max_val) ? max_val : val_i);
      val_j = (val_j < min_val) ? min_val : ((val_j > max_val) ? max_val : val_j);
      
      (*out_mem) = val_i;  out_mem++;
      (*out_mem) = val_j;  out_mem++;
      }
    
    const uword i = j-1;
    
    if(i < N)
      {
      eT val_i = A[i];
      
      val_i = (val_i < min_val) ? min_val : ((val_i > max_val) ? max_val : val_i);
      
      (*out_mem) = val_i;
      }
    }
  else
    {
    for(uword col=0; col<n_cols; ++col)
    for(uword row=0; row<n_rows; ++row)
      {
      eT val = P.at(row,col);
      
      val = (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
      
      (*out_mem) = val;  out_mem++;
      }
    }
  }



//! @}
