// Copyright (C) 2010-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_cross
//! @{



//! cross product (only valid for 3 dimensional vectors)
template<typename T1, typename T2>
arma_warn_unused
inline
const Glue<T1, T2, glue_cross>
cross(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_cross>(X.get_ref(), Y.get_ref());
  }



//! @}
