// Copyright (C) 2016 National ICT Australia (NICTA)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
//
// Written by Yixuan Qiu


namespace newarp
{


//! Tiny functions to check attributes of complex numbers
struct cx_attrib
  {
  template<typename T>
  arma_inline static bool is_real   (const std::complex<T>& v, const T eps) { return (std::abs(v.imag()) <= eps); }
  
  template<typename T>
  arma_inline static bool is_complex(const std::complex<T>& v, const T eps) { return (std::abs(v.imag()) >  eps); }
  
  template<typename T>
  arma_inline static bool is_conj(const std::complex<T>& v1, const std::complex<T>& v2, const T eps)  { return (std::abs(v1 - std::conj(v2)) <= eps); }
  };


}  // namespace newarp
