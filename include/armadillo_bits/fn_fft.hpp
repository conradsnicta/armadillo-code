// Copyright (C) 2013-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_fft
//! @{



// 1D FFT & 1D IFFT



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_real<typename T1::elem_type>::value),
  const mtOp<std::complex<typename T1::pod_type>, T1, op_fft_real>
  >::result
fft(const T1& A)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<std::complex<typename T1::pod_type>, T1, op_fft_real>(A, uword(0), uword(1));
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_real<typename T1::elem_type>::value),
  const mtOp<std::complex<typename T1::pod_type>, T1, op_fft_real>
  >::result
fft(const T1& A, const uword N)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<std::complex<typename T1::pod_type>, T1, op_fft_real>(A, N, uword(0));
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_complex_strict<typename T1::elem_type>::value),
  const Op<T1, op_fft_cx>
  >::result
fft(const T1& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_fft_cx>(A, uword(0), uword(1));
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_complex_strict<typename T1::elem_type>::value),
  const Op<T1, op_fft_cx>
  >::result
fft(const T1& A, const uword N)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_fft_cx>(A, N, uword(0));
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_complex_strict<typename T1::elem_type>::value),
  const Op<T1, op_ifft_cx>
  >::result
ifft(const T1& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_ifft_cx>(A, uword(0), uword(1));
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_complex_strict<typename T1::elem_type>::value),
  const Op<T1, op_ifft_cx>
  >::result
ifft(const T1& A, const uword N)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_ifft_cx>(A, N, uword(0));
  }



//! @}
