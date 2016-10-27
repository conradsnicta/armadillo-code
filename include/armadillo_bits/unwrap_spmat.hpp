// Copyright (C) 2012-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup unwrap_spmat
//! @{



template<typename T1>
struct unwrap_spmat
  {
  typedef typename T1::elem_type eT;
  
  typedef SpMat<eT> stored_type;
  
  inline
  unwrap_spmat(const T1& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<eT> M;
  };



template<typename eT>
struct unwrap_spmat< SpMat<eT> >
  {
  typedef SpMat<eT> stored_type;
  
  inline
  unwrap_spmat(const SpMat<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<eT>& M;
  };



template<typename eT>
struct unwrap_spmat< SpRow<eT> >
  {
  typedef SpRow<eT> stored_type;
  
  inline
  unwrap_spmat(const SpRow<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpRow<eT>& M;
  };



template<typename eT>
struct unwrap_spmat< SpCol<eT> >
  {
  typedef SpCol<eT> stored_type;
  
  inline
  unwrap_spmat(const SpCol<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpCol<eT>& M;
  };



template<typename T1, typename spop_type>
struct unwrap_spmat< SpOp<T1, spop_type> >
  {
  typedef typename T1::elem_type eT;
  
  typedef SpMat<eT> stored_type;
  
  inline
  unwrap_spmat(const SpOp<T1, spop_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<eT> M;
  };



template<typename T1, typename T2, typename spglue_type>
struct unwrap_spmat< SpGlue<T1, T2, spglue_type> >
  {
  typedef typename T1::elem_type eT;
  
  typedef SpMat<eT> stored_type;
  
  inline
  unwrap_spmat(const SpGlue<T1, T2, spglue_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<eT> M;
  };



template<typename out_eT, typename T1, typename spop_type>
struct unwrap_spmat< mtSpOp<out_eT, T1, spop_type> >
  {
  typedef SpMat<out_eT> stored_type;
  
  inline
  unwrap_spmat(const mtSpOp<out_eT, T1, spop_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<out_eT> M;
  };



//! @}
