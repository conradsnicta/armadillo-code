// Copyright (C) 2008-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_eye
//! @{



arma_warn_unused
arma_inline
const Gen<mat, gen_eye>
eye(const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  return Gen<mat, gen_eye>(n_rows, n_cols);
  }



arma_warn_unused
arma_inline
const Gen<mat, gen_eye>
eye(const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  return Gen<mat, gen_eye>(s.n_rows, s.n_cols);
  }



template<typename obj_type>
arma_warn_unused
arma_inline
const Gen<obj_type, gen_eye>
eye(const uword n_rows, const uword n_cols, const typename arma_Mat_Col_Row_only<obj_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(is_Col<obj_type>::value)
    {
    arma_debug_check( (n_cols != 1), "eye(): incompatible size" );
    }
  else
  if(is_Row<obj_type>::value)
    {
    arma_debug_check( (n_rows != 1), "eye(): incompatible size" );
    }
  
  return Gen<obj_type, gen_eye>(n_rows, n_cols);
  }



template<typename obj_type>
arma_warn_unused
arma_inline
const Gen<obj_type, gen_eye>
eye(const SizeMat& s, const typename arma_Mat_Col_Row_only<obj_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return eye<obj_type>(s.n_rows, s.n_cols);
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
eye(const uword n_rows, const uword n_cols, const typename arma_SpMat_SpCol_SpRow_only<obj_type>::result* junk = NULL)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(is_SpCol<obj_type>::value == true)
    {
    arma_debug_check( (n_cols != 1), "eye(): incompatible size" );
    }
  else
  if(is_SpRow<obj_type>::value == true)
    {
    arma_debug_check( (n_rows != 1), "eye(): incompatible size" );
    }
  
  obj_type out;
  
  out.eye(n_rows, n_cols);
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
eye(const SizeMat& s, const typename arma_SpMat_SpCol_SpRow_only<obj_type>::result* junk = NULL)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return eye<obj_type>(s.n_rows, s.n_cols);
  }



//! @}
