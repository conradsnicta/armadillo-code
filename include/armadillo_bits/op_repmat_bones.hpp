// Copyright (C) 2009-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by Dimitrios Bouzas



//! \addtogroup op_repmat
//! @{



class op_repmat
  {
  public:
  
  template<typename obj> inline static void apply_noalias(Mat<typename obj::elem_type>& out, const obj& X, const uword copies_per_row, const uword copies_per_col);
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_repmat>& in);
  };



//! @}
