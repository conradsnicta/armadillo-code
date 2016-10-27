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


TEST_CASE("attributes_1")
  {
  mat A(5,6);
  REQUIRE(A.n_rows == 5);
  REQUIRE(A.n_cols == 6);
  REQUIRE(A.n_elem == 30);
  
  vec B(5);
  REQUIRE(B.n_rows == 5);
  REQUIRE(B.n_cols == 1);
  REQUIRE(B.n_elem == 5);
  
  rowvec C(6);
  REQUIRE(C.n_rows == 1);
  REQUIRE(C.n_cols == 6);
  REQUIRE(C.n_elem == 6);
  
  cube D(5,6,2);
  REQUIRE(D.n_rows   == 5);
  REQUIRE(D.n_cols   == 6);
  REQUIRE(D.n_slices == 2);
  REQUIRE(D.n_elem   == 60);
  
  sp_mat E(50,60);
  E(0,0) = 1.0;
  E(E.n_rows-1,E.n_cols-1) = 1.0;
  
  REQUIRE(E.n_rows    == 50);
  REQUIRE(E.n_cols    == 60);
  REQUIRE(E.n_elem    == 3000);
  REQUIRE(E.n_nonzero == 2);
  
  field<mat> G(5,6,2);
  REQUIRE(G.n_rows   == 5);
  REQUIRE(G.n_cols   == 6);
  REQUIRE(G.n_slices == 2);
  REQUIRE(G.n_elem   == 60);
  }



