// Copyright (C) 2012-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_spones
//! @{



//! Generate a sparse matrix with the non-zero values in the same locations as in the given sparse matrix X,
//! with the non-zero values set to one
template<typename T1>
arma_warn_unused
inline
SpMat<typename T1::elem_type>
spones(const SpBase<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  SpMat<eT> out( X.get_ref() );
  
  arrayops::inplace_set( access::rwp(out.values), eT(1), out.n_nonzero );
  
  return out;
  }



//! @}
