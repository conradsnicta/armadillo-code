// Copyright (C) 2011-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup arrayops
//! @{


class arrayops
  {
  public:
  
  template<typename eT>
  arma_hot arma_inline static void
  copy(eT* dest, const eT* src, const uword n_elem);
  
  
  template<typename eT>
  arma_hot inline static void
  copy_small(eT* dest, const eT* src, const uword n_elem);
  
  
  template<typename eT>
  arma_hot inline static void
  copy_forwards(eT* dest, const eT* src, const uword n_elem);
  
  
  template<typename eT>
  arma_hot inline static void
  copy_backwards(eT* dest, const eT* src, const uword n_elem);
  
  
  template<typename eT>
  arma_hot inline static void
  fill_zeros(eT* dest, const uword n_elem);
  
  
  template<typename eT>
  arma_hot inline static void
  replace(eT* mem, const uword n_elem, const eT old_val, const eT new_val);
  
  
  // 
  // array = convert(array)
  
  template<typename out_eT, typename in_eT>
  arma_hot arma_inline static void
  convert_cx_scalar(out_eT& out, const in_eT&  in, const typename arma_not_cx<out_eT>::result* junk1 = 0, const typename arma_not_cx< in_eT>::result* junk2 = 0);
  
  template<typename out_eT, typename in_T>
  arma_hot arma_inline static void
  convert_cx_scalar(out_eT& out, const std::complex<in_T>& in, const typename arma_not_cx<out_eT>::result* junk = 0);
  
  template<typename out_T, typename in_T>
  arma_hot arma_inline static void
  convert_cx_scalar(std::complex<out_T>& out, const std::complex< in_T>& in);
  
  template<typename out_eT, typename in_eT>
  arma_hot inline static void
  convert(out_eT* dest, const in_eT* src, const uword n_elem);
  
  template<typename out_eT, typename in_eT>
  arma_hot inline static void
  convert_cx(out_eT* dest, const in_eT* src, const uword n_elem);
  
  
  // 
  // array op= array
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_plus(eT* dest, const eT* src, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_minus(eT* dest, const eT* src, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_mul(eT* dest, const eT* src, const uword n_elem);
   
  template<typename eT>
  arma_hot inline static
  void
  inplace_div(eT* dest, const eT* src, const uword n_elem);
  
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_plus_base(eT* dest, const eT* src, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_minus_base(eT* dest, const eT* src, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_mul_base(eT* dest, const eT* src, const uword n_elem);
   
  template<typename eT>
  arma_hot inline static
  void
  inplace_div_base(eT* dest, const eT* src, const uword n_elem);
  
  
  // 
  // array op= scalar
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_set(eT* dest, const eT val, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_set_base(eT* dest, const eT val, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_set_small(eT* dest, const eT val, const uword n_elem);
  
  template<typename eT, const uword n_elem>
  arma_hot inline static
  void
  inplace_set_fixed(eT* dest, const eT val);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_plus(eT* dest, const eT val, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_minus(eT* dest, const eT val, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static void
  inplace_mul(eT* dest, const eT val, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_div(eT* dest, const eT val, const uword n_elem);
  
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_plus_base(eT* dest, const eT val, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_minus_base(eT* dest, const eT val, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static void
  inplace_mul_base(eT* dest, const eT val, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_div_base(eT* dest, const eT val, const uword n_elem);
  
  
  // 
  // scalar = op(array)
  
  template<typename eT>
  arma_hot inline static
  eT
  accumulate(const eT* src, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  eT
  product(const eT* src, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  bool
  is_finite(const eT* src, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  bool
  has_inf(const eT* src, const uword n_elem);
  
  template<typename eT>
  arma_hot inline static
  bool
  has_nan(const eT* src, const uword n_elem);
  };



//! @}
