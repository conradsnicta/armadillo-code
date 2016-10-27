// Copyright (C) 2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


#include <armadillo>
#include "catch.hpp"

using namespace arma;


TEST_CASE("instantiation_mat_1")
  {
  const uword n_rows = 5;
  const uword n_cols = 6;
  
   mat A(n_rows,n_cols);
  fmat B(n_rows,n_cols);
  umat C(n_rows,n_cols);
  imat D(n_rows,n_cols);
  
  cx_mat  E(n_rows,n_cols);
  cx_fmat F(n_rows,n_cols);
  
  mat::fixed<5,6> G;
  }


// TODO: colvec_instantiation
// TODO: rowvec_instantiation


TEST_CASE("instantiation_cube_1")
  {
  const uword n_rows   = 5;
  const uword n_cols   = 6;
  const uword n_slices = 2;
  
   cube A(n_rows,n_cols,n_slices);
  fcube B(n_rows,n_cols,n_slices);
  ucube C(n_rows,n_cols,n_slices);
  icube D(n_rows,n_cols,n_slices);
  
  cx_cube  E(n_rows,n_cols,n_slices);
  cx_fcube F(n_rows,n_cols,n_slices);
  
  cube::fixed<5,6,2> G;
  }


// TODO: field_instantiation 1D 2D 3D
// TODO: spmat_instantiation


