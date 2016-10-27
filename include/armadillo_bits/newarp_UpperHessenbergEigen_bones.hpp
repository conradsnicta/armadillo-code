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


//! Calculate the eigenvalues and eigenvectors of an upper Hessenberg matrix.
//! This class is a wrapper of the Lapack functions `_lahqr` and `_trevc`.
template<typename eT>
class UpperHessenbergEigen
  {
  private:

  blas_int                n;
  Mat<eT>                 mat_Z;  // In the first stage, H = ZTZ', Z is an orthogonal matrix
                                  // In the second stage, Z will be overwritten by the eigenvectors of H
  Mat<eT>                 mat_T;  // H = ZTZ', T is a Schur form matrix
  Col< std::complex<eT> > evals;  // eigenvalues of H
  bool                    computed;


  public:

  //! Default constructor. Computation can
  //! be performed later by calling the compute() method.
  inline UpperHessenbergEigen();

  //! Constructor to create an object that calculates the eigenvalues
  //! and eigenvectors of an upper Hessenberg matrix `mat_obj`.
  inline UpperHessenbergEigen(const Mat<eT>& mat_obj);

  //! Compute the eigenvalue decomposition of an upper Hessenberg matrix.
  inline void compute(const Mat<eT>& mat_obj);

  //! Retrieve the eigenvalues.
  inline Col< std::complex<eT> > eigenvalues();

  //! Retrieve the eigenvectors.
  inline Mat< std::complex<eT> > eigenvectors();
  };


}  // namespace newarp
