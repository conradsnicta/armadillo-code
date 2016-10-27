// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_nonzeros
//! @{


template<typename T1>
arma_warn_unused
inline
const Op<T1,op_nonzeros>
nonzeros(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1,op_nonzeros>(X.get_ref());
  }



template<typename T1>
arma_warn_unused
inline
Col<typename T1::elem_type>
nonzeros(const SpBase<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Col<typename T1::elem_type> out;
  
  op_nonzeros::apply_noalias(out, X.get_ref());
  
  return out;
  }



//! @}
