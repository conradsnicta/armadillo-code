// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_eig_gen
//! @{


template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, Col< std::complex<typename T1::pod_type> > >::result
eig_gen
  (
  const Base<typename T1::elem_type, T1>& expr
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type     T;
  typedef typename std::complex<T> eT;
  
  Col<eT> eigvals;
  Mat<eT> eigvecs;
  
  const bool status = auxlib::eig_gen(eigvals, eigvecs, false, expr.get_ref());
  
  if(status == false)
    {
    eigvals.reset();
    arma_stop_runtime_error("eig_gen(): decomposition failed");
    }
  
  return eigvals;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, bool >::result
eig_gen
  (
         Col< std::complex<typename T1::pod_type> >& eigvals,
  const Base< typename T1::elem_type, T1>&           expr
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type     T;
  typedef typename std::complex<T> eT;
  
  Mat<eT> eigvecs;
  
  const bool status = auxlib::eig_gen(eigvals, eigvecs, false, expr.get_ref());
  
  if(status == false)
    {
    eigvals.reset();
    arma_debug_warn("eig_gen(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, bool >::result
eig_gen
  (
        Col< std::complex<typename T1::pod_type> >& eigvals,
        Mat< std::complex<typename T1::pod_type> >& eigvecs,
  const Base<typename T1::elem_type, T1>&           expr
  )
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (void_ptr(&eigvals) == void_ptr(&eigvecs)), "eig_gen(): parameter 'eigval' is an alias of parameter 'eigvec'" );
  
  const bool status = auxlib::eig_gen(eigvals, eigvecs, true, expr.get_ref());
  
  if(status == false)
    {
    eigvals.reset();
    eigvecs.reset();
    arma_debug_warn("eig_gen(): decomposition failed");
    }
  
  return status;
  }



//! @}
