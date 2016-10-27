// Copyright (C) 2012-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup op_hist
//! @{



class op_hist
  {
  public:
  
  template<typename eT>
  inline static void apply_noalias(Mat<uword>& out, const Mat<eT>& A, const uword n_bins, const bool A_is_row);
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const mtOp<uword, T1, op_hist>& X);
  };



//! @}
