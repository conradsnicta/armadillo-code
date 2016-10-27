// Copyright (C) 2008-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_sort_index
//! @{



class op_sort_index
  {
  public:
  
  template<typename T1>
  static inline bool apply_noalias(Mat<uword>& out, const Proxy<T1>& P, const uword sort_type);
  
  template<typename T1>
  static inline void apply(Mat<uword>& out, const mtOp<uword,T1,op_sort_index>& in);
  };



class op_stable_sort_index
  {
  public:
  
  template<typename T1>
  static inline bool apply_noalias(Mat<uword>& out, const Proxy<T1>& P, const uword sort_type);
  
  template<typename T1>
  static inline void apply(Mat<uword>& out, const mtOp<uword,T1,op_stable_sort_index>& in);
  };



template<typename eT>
struct arma_sort_index_packet
  {
  eT    val;
  uword index;
  };



template<typename eT>
struct arma_sort_index_helper_ascend
  {
  arma_inline
  bool
  operator() (const arma_sort_index_packet<eT>& A, const arma_sort_index_packet<eT>& B) const
    {
    return (A.val < B.val);
    }
  };



template<typename eT>
struct arma_sort_index_helper_descend
  {
  arma_inline
  bool
  operator() (const arma_sort_index_packet<eT>& A, const arma_sort_index_packet<eT>& B) const
    {
    return (A.val > B.val);
    }
  };



template<typename T>
struct arma_sort_index_helper_ascend< std::complex<T> >
  {
  typedef typename std::complex<T> eT;
  
  inline
  bool
  operator() (const arma_sort_index_packet<eT>& A, const arma_sort_index_packet<eT>& B) const
    {
    return (std::abs(A.val) < std::abs(B.val));
    }
  
  // inline
  // bool
  // operator() (const arma_sort_index_packet<eT>& A, const arma_sort_index_packet<eT>& B) const
  //   {
  //   const T abs_A_val = std::abs(A.val);
  //   const T abs_B_val = std::abs(B.val);
  //   
  //   return ( (abs_A_val != abs_B_val) ? (abs_A_val < abs_B_val) : (std::arg(A.val) < std::arg(B.val)) );
  //   }
  };



template<typename T>
struct arma_sort_index_helper_descend< std::complex<T> >
  {
  typedef typename std::complex<T> eT;
  
  inline
  bool
  operator() (const arma_sort_index_packet<eT>& A, const arma_sort_index_packet<eT>& B) const
    {
    return (std::abs(A.val) > std::abs(B.val));
    }
  
  // inline
  // bool
  // operator() (const arma_sort_index_packet<eT>& A, const arma_sort_index_packet<eT>& B) const
  //   {
  //   const T abs_A_val = std::abs(A.val);
  //   const T abs_B_val = std::abs(B.val);
  //   
  //   return ( (abs_A_val != abs_B_val) ? (abs_A_val > abs_B_val) : (std::arg(A.val) > std::arg(B.val)) );
  //   }
  };



//! @}
