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


//! This class implements the eigen solver for real symmetric matrices.
template<typename eT, int SelectionRule, typename OpType>
class SymEigsSolver
  {
  protected:

  const OpType&     op;        // object to conduct matrix operation, e.g. matrix-vector product
  const uword       nev;       // number of eigenvalues requested
  Col<eT>           ritz_val;  // ritz values

  // Sort the first nev Ritz pairs in decreasing magnitude order
  // This is used to return the final results
  virtual void sort_ritzpair();


  private:

  const uword       dim_n;     // dimension of matrix A
  const uword       ncv;       // number of ritz values
  uword             nmatop;    // number of matrix operations called
  uword             niter;     // number of restarting iterations
  Mat<eT>           fac_V;     // V matrix in the Arnoldi factorisation
  Mat<eT>           fac_H;     // H matrix in the Arnoldi factorisation
  Col<eT>           fac_f;     // residual in the Arnoldi factorisation
  Mat<eT>           ritz_vec;  // ritz vectors
  std::vector<bool> ritz_conv; // indicator of the convergence of ritz values
  const eT          prec;      // precision parameter used to test convergence
                               // prec = epsilon^(2/3)
                               // epsilon is the machine precision,
                               // e.g. ~= 1e-16 for the "double" type

  // Arnoldi factorisation starting from step-k
  inline void factorise_from(uword from_k, uword to_m, const Col<eT>& fk);

  // Implicitly restarted Arnoldi factorisation
  inline void restart(uword k);

  // Calculate the number of converged Ritz values
  inline uword num_converged(eT tol);

  // Return the adjusted nev for restarting
  inline uword nev_adjusted(uword nconv);

  // Retrieve and sort ritz values and ritz vectors
  inline void retrieve_ritzpair();


  public:

  //! Constructor to create a solver object.
  inline SymEigsSolver(const OpType& op_, uword nev_, uword ncv_);

  //! Providing the initial residual vector for the algorithm.
  inline void init(eT* init_resid);

  //! Providing a random initial residual vector.
  inline void init();

  //! Conducting the major computation procedure.
  inline uword compute(uword maxit = 1000, eT tol = 1e-10);

  //! Returning the number of iterations used in the computation.
  inline uword num_iterations() { return niter; }

  //! Returning the number of matrix operations used in the computation.
  inline uword num_operations() { return nmatop; }

  //! Returning the converged eigenvalues.
  inline Col<eT> eigenvalues();

  //! Returning the eigenvectors associated with the converged eigenvalues.
  inline Mat<eT> eigenvectors(uword nvec);
  //! Returning all converged eigenvectors.
  inline Mat<eT> eigenvectors() { return eigenvectors(nev); }
  };


}  // namespace newarp
