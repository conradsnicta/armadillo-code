// Copyright (C) 2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_find_unique
//! @{



class op_find_unique
  {
  public:
  
  template<typename T1>
  static inline bool apply_helper(Mat<uword>& out, const Proxy<T1>& P, const bool ascending_indices);
  
  template<typename T1>
  static inline void apply(Mat<uword>& out, const mtOp<uword,T1,op_find_unique>& in);
  };



template<typename eT>
struct arma_find_unique_packet
  {
  eT    val;
  uword index;
  };



template<typename eT>
struct arma_find_unique_comparator
  {
  arma_inline
  bool
  operator() (const arma_find_unique_packet<eT>& A, const arma_find_unique_packet<eT>& B) const
    {
    return (A.val < B.val);
    }
  };



template<typename T>
struct arma_find_unique_comparator< std::complex<T> >
  {
  arma_inline
  bool
  operator() (const arma_find_unique_packet< std::complex<T> >& A, const arma_find_unique_packet< std::complex<T> >& B) const
    {
    const T A_real = A.val.real();
    const T B_real = B.val.real();
    
    return (  (A_real < B_real) ? true : ((A_real == B_real) ? (A.val.imag() < B.val.imag()) : false)  );
    }
  };



//! @}
