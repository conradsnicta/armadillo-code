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


template<typename eT>
inline 
DenseGenMatProd<eT>::DenseGenMatProd(const Mat<eT>& mat_obj)
  : op_mat(mat_obj)
  , n_rows(mat_obj.n_rows)
  , n_cols(mat_obj.n_cols)
  {
  arma_extra_debug_sigprint();
  }



// Perform the matrix-vector multiplication operation \f$y=Ax\f$.
// y_out = A * x_in
template<typename eT>
inline
void
DenseGenMatProd<eT>::perform_op(eT* x_in, eT* y_out) const
  {
  arma_extra_debug_sigprint();
  
  Col<eT> x(x_in , n_cols, false);
  Col<eT> y(y_out, n_rows, false);
  y = op_mat * x;
  }


}  // namespace newarp
