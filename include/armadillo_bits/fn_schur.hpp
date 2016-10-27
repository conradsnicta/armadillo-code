// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_schur
//! @{


template<typename T1>
inline
bool
schur
  (
         Mat<typename T1::elem_type>&    S,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> U;
  
  const bool status = auxlib::schur(U, S, X.get_ref(), false);
  
  if(status == false)
    {
    S.reset();
    arma_debug_warn("schur(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
arma_warn_unused
inline
Mat<typename T1::elem_type>
schur
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> S;
  Mat<eT> U;
  
  const bool status = auxlib::schur(U, S, X.get_ref(), false);
  
  if(status == false)
    {
    S.reset();
    arma_stop_runtime_error("schur(): decomposition failed");
    }
  
  return S;
  }



template<typename T1> 
inline
bool
schur
  (
         Mat<typename T1::elem_type>&    U,
         Mat<typename T1::elem_type>&    S,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check( void_ptr(&U) == void_ptr(&S), "schur(): 'U' is an alias of 'S'" );
  
  const bool status = auxlib::schur(U, S, X.get_ref(), true);
  
  if(status == false)
    {
    U.reset();
    S.reset();
    arma_debug_warn("schur(): decomposition failed");
    }
  
  return status;
  }



//! @}
