// Copyright (C) 2009-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_prod
//! @{


class op_prod
  {
  public:
  
  template<typename eT>
  inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1, op_prod>& in);
  
  template<typename eT>
  inline static eT prod(const subview<eT>& S);
  
  template<typename T1>
  inline static typename T1::elem_type prod(const Base<typename T1::elem_type,T1>& X);
  };


//! @}
