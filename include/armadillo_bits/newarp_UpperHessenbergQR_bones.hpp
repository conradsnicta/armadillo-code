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


//! Perform the QR decomposition of an upper Hessenberg matrix.
template<typename eT>
class UpperHessenbergQR
  {
  protected:

  uword   n;
  Mat<eT> mat_T;
  // Gi = [ cos[i]  sin[i]]
  //      [-sin[i]  cos[i]]
  // Q = G1 * G2 * ... * G_{n-1}
  Col<eT> rot_cos;
  Col<eT> rot_sin;
  bool    computed;


  public:

  //! Default constructor. Computation can
  //! be performed later by calling the compute() method.
  inline UpperHessenbergQR();

  //! Constructor to create an object that performs and stores the
  //! QR decomposition of an upper Hessenberg matrix `mat_obj`.
  inline UpperHessenbergQR(const Mat<eT>& mat_obj);

  //! Conduct the QR factorisation of an upper Hessenberg matrix.
  virtual void compute(const Mat<eT>& mat_obj);

  //! Return the \f$RQ\f$ matrix, the multiplication of \f$R\f$ and \f$Q\f$,
  //! which is an upper Hessenberg matrix.
  virtual Mat<eT> matrix_RQ();

  //! Apply the \f$Q\f$ matrix to another matrix \f$Y\f$.
  inline void apply_YQ(Mat<eT>& Y);
  };



//! Perform the QR decomposition of a tridiagonal matrix, a special
//! case of upper Hessenberg matrices.
template<typename eT>
class TridiagQR : public UpperHessenbergQR<eT>
  {
  public:

  //! Default constructor. Computation can
  //! be performed later by calling the compute() method.
  inline TridiagQR();

  //! Constructor to create an object that performs and stores the
  //! QR decomposition of a tridiagonal matrix `mat_obj`.
  inline TridiagQR(const Mat<eT>& mat_obj);

  //! Conduct the QR factorisation of a tridiagonal matrix.
  inline void compute(const Mat<eT>& mat_obj);

  //! Return the \f$RQ\f$ matrix, the multiplication of \f$R\f$ and \f$Q\f$,
  //! which is a tridiagonal matrix.
  inline Mat<eT> matrix_RQ();
  };


}  // namespace newarp
