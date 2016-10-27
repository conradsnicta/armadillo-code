// Copyright (C) 2009-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by Dimitrios Bouzas
// Written by Stanislav Funiak


//! \addtogroup fn_rank
//! @{



template<typename T1>
arma_warn_unused
inline
uword
rank
  (
  const Base<typename T1::elem_type,T1>& X,
        typename T1::pod_type            tol = 0.0,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  uword  X_n_rows;
  uword  X_n_cols;
  Col<T> s;
  
  const bool status = auxlib::svd_dc(s, X, X_n_rows, X_n_cols);
  
  if(status == false)
    {
    arma_stop_runtime_error("rank(): svd failed");
    
    return uword(0);
    }
  
  const uword s_n_elem = s.n_elem;
  const T*    s_mem    = s.memptr();
  
  // set tolerance to default if it hasn't been specified
  if( (tol == T(0)) && (s_n_elem > 0) )
    {
    tol = (std::max)(X_n_rows, X_n_cols) * s_mem[0] * std::numeric_limits<T>::epsilon();
    }
  
  uword count = 0;
  
  for(uword i=0; i < s_n_elem; ++i)  { count += (s_mem[i] > tol) ? uword(1) : uword(0); }
  
  return count;
  }



//! @}
