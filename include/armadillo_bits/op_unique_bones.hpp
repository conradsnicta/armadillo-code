// Copyright (C) 2012-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by Arnold Wiliem



//! \addtogroup op_unique
//! @{



class op_unique
  {
  public:
  
  template<typename T1>
  inline static bool apply_helper(Mat<typename T1::elem_type>& out, const Proxy<T1>& P);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_unique>& in);
  };



template<typename eT>
struct arma_unique_comparator
  {
  arma_inline
  bool
  operator() (const eT a, const eT b) const
    {
    return ( a < b );
    }
  };



template<typename T>
struct arma_unique_comparator< std::complex<T> >
  {
  arma_inline
  bool
  operator() (const std::complex<T>& a, const std::complex<T>& b) const
    {
    const T a_real = a.real();
    const T b_real = b.real();
    
    return (  (a_real < b_real) ? true : ((a_real == b_real) ? (a.imag() < b.imag()) : false)  );
    }
  };



//! @}
