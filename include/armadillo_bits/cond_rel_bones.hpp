// Copyright (C) 2012 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup cond_rel
//! @{


//
// for preventing pedantic compiler warnings

template<const bool do_eval>
class cond_rel
  {
  public:
  
  template<typename eT> arma_inline static bool lt(const eT A, const eT B);
  template<typename eT> arma_inline static bool gt(const eT A, const eT B);

  template<typename eT> arma_inline static bool leq(const eT A, const eT B);
  template<typename eT> arma_inline static bool geq(const eT A, const eT B);
  
  template<typename eT> arma_inline static eT make_neg(const eT val);
  };



//! @}
