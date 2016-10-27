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
class DoubleShiftQR
  {
  private:

  uword               n;        // Dimension of the matrix
  Mat<eT>             mat_H;    // A copy of the matrix to be factorised
  eT                  shift_s;  // Shift constant
  eT                  shift_t;  // Shift constant
  Mat<eT>             ref_u;    // Householder reflectors
  Col<unsigned short> ref_nr;   // How many rows does each reflector affects
                                // 3 - A general reflector
                                // 2 - A Givens rotation
                                // 1 - An identity transformation
  const eT            prec;     // Approximately zero
  const eT            eps_rel;
  const eT            eps_abs;
  bool                computed; // Whether matrix has been factorised

  inline      void compute_reflector(const eT& x1, const eT& x2, const eT& x3, uword ind);
  arma_inline void compute_reflector(const eT* x, uword ind);

  // Update the block X = H(il:iu, il:iu)
  inline void update_block(uword il, uword iu);

  // P = I - 2 * u * u' = P'
  // PX = X - 2 * u * (u'X)
  inline void apply_PX(Mat<eT>& X, uword oi, uword oj, uword nrow, uword ncol, uword u_ind);

  // x is a pointer to a vector
  // Px = x - 2 * dot(x, u) * u
  inline void apply_PX(eT* x, uword u_ind);

  // XP = X - 2 * (X * u) * u'
  inline void apply_XP(Mat<eT>& X, uword oi, uword oj, uword nrow, uword ncol, uword u_ind);


  public:

  inline DoubleShiftQR(uword size);

  inline DoubleShiftQR(const Mat<eT>& mat_obj, eT s, eT t);

  inline void compute(const Mat<eT>& mat_obj, eT s, eT t);

  inline Mat<eT> matrix_QtHQ();

  inline void apply_QtY(Col<eT>& y);

  inline void apply_YQ(Mat<eT>& Y);
  };


}  // namespace newarp
