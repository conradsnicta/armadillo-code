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


//! Define matrix operations on existing matrix objects
template<typename eT>
class SparseGenMatProd
  {
  private:
  
  const SpMat<eT>& op_mat;
  
  
  public:
  
  const uword n_rows;  // number of rows of the underlying matrix
  const uword n_cols;  // number of columns of the underlying matrix
  
  inline SparseGenMatProd(const SpMat<eT>& mat_obj);
  
  inline void perform_op(eT* x_in, eT* y_out) const;
  };


}  // namespace newarp
