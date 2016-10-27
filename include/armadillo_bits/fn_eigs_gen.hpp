// Copyright (C) 2013-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by Ryan Curtin


//! \addtogroup fn_eigs_gen
//! @{


//! eigenvalues of general sparse matrix X
template<typename T1>
arma_warn_unused
inline
Col< std::complex<typename T1::pod_type> >
eigs_gen
  (
  const SpBase<typename T1::elem_type, T1>& X,
  const uword                               n_eigvals,
  const char*                               form = "lm",
  const typename T1::pod_type               tol  = 0.0,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  Mat< std::complex<T> > eigvec;
  Col< std::complex<T> > eigval;
  
  const bool status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, form, tol);
  
  if(status == false)
    {
    eigval.reset();
    arma_stop_runtime_error("eigs_gen(): decomposition failed");
    }
  
  return eigval;
  }



//! eigenvalues of general sparse matrix X
template<typename T1>
inline
bool
eigs_gen
  (
           Col< std::complex<typename T1::pod_type> >& eigval,
  const SpBase<typename T1::elem_type, T1>&            X,
  const uword                                          n_eigvals,
  const char*                                          form = "lm",
  const typename T1::pod_type                          tol  = 0.0,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  Mat< std::complex<T> > eigvec;
  
  const bool status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, form, tol);
  
  if(status == false)
    {
    eigval.reset();
    arma_debug_warn("eigs_gen(): decomposition failed");
    }
  
  return status;
  }



//! eigenvalues and eigenvectors of general real sparse matrix X
template<typename T1>
inline
bool
eigs_gen
  (
         Col< std::complex<typename T1::pod_type> >& eigval,
         Mat< std::complex<typename T1::pod_type> >& eigvec,
  const SpBase<typename T1::elem_type, T1>&          X,
  const uword                                        n_eigvals,
  const char*                                        form = "lm",
  const typename T1::pod_type                        tol  = 0.0,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check( void_ptr(&eigval) == void_ptr(&eigvec), "eigs_gen(): parameter 'eigval' is an alias of parameter 'eigvec'" );
  
  const bool status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, form, tol);
  
  if(status == false)
    {
    eigval.reset();
    eigvec.reset();
    arma_debug_warn("eigs_gen(): decomposition failed");
    }
  
  return status;
  }



//! @}
