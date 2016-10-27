// Copyright (C) 2012-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by Ryan Curtin


//! \addtogroup fn_speye
//! @{



//! Generate a sparse matrix with the values along the main diagonal set to one
template<typename obj_type>
arma_warn_unused
inline
obj_type
speye(const uword n_rows, const uword n_cols, const typename arma_SpMat_SpCol_SpRow_only<obj_type>::result* junk = NULL)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(is_SpCol<obj_type>::value == true)
    {
    arma_debug_check( (n_cols != 1), "speye(): incompatible size" );
    }
  else
  if(is_SpRow<obj_type>::value == true)
    {
    arma_debug_check( (n_rows != 1), "speye(): incompatible size" );
    }
  
  obj_type out;
  
  out.eye(n_rows, n_cols);
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
speye(const SizeMat& s, const typename arma_SpMat_SpCol_SpRow_only<obj_type>::result* junk = NULL)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return speye<obj_type>(s.n_rows, s.n_cols);
  }



// Convenience shortcut method (no template parameter necessary)
arma_warn_unused
inline
sp_mat
speye(const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  sp_mat out;
  
  out.eye(n_rows, n_cols);
  
  return out;
  }



arma_warn_unused
inline
sp_mat
speye(const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  sp_mat out;
  
  out.eye(s.n_rows, s.n_cols);
  
  return out;
  }



//! @}
