// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup GenSpecialiser
//! @{


template<typename elem_type, bool is_gen_zeros, bool is_gen_ones, bool is_gen_randu, bool is_gen_randn>
struct GenSpecialiser
  {
  arma_inline elem_type generate() const { return elem_type(); }
  };


template<typename elem_type>
struct GenSpecialiser<elem_type, true, false, false, false>
  {
  arma_inline elem_type generate() const { return elem_type(0); }
  };


template<typename elem_type>
struct GenSpecialiser<elem_type, false, true, false, false>
  {
  arma_inline elem_type generate() const { return elem_type(1); }
  };


template<typename elem_type>
struct GenSpecialiser<elem_type, false, false, true, false>
  {
  arma_inline elem_type generate() const { return elem_type(arma_rng::randu<elem_type>()); }
  };


template<typename elem_type>
struct GenSpecialiser<elem_type, false, false, false, true>
  {
  arma_inline elem_type generate() const { return elem_type(arma_rng::randn<elem_type>()); }
  };


//! @}
