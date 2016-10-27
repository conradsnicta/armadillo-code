// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_trapz
//! @{



template<typename T1, typename T2>
arma_warn_unused
inline
const Glue<T1, T2, glue_trapz>
trapz
  (
  const Base<typename T1::elem_type,T1>& X,
  const Base<typename T1::elem_type,T2>& Y,
  const uword                            dim = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_trapz>(X.get_ref(), Y.get_ref(), dim);
  }



template<typename T1>
arma_warn_unused
inline
const Op<T1, op_trapz>
trapz
  (
  const Base<typename T1::elem_type,T1>& Y,
  const uword                            dim = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_trapz>(Y.get_ref(), dim, uword(0));
  }



//! @}
