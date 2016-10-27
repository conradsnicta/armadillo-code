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


TEST_CASE("gen_linspace_1")
  {
  vec a = linspace(1,5,5);
  
  REQUIRE(a(0) == Approx(1.0));
  REQUIRE(a(1) == Approx(2.0));
  REQUIRE(a(2) == Approx(3.0));
  REQUIRE(a(3) == Approx(4.0));
  REQUIRE(a(4) == Approx(5.0));
  
  vec b = linspace<vec>(1,5,6);
  
  REQUIRE(b(0) == Approx(1.0));
  REQUIRE(b(1) == Approx(1.8));
  REQUIRE(b(2) == Approx(2.6));
  REQUIRE(b(3) == Approx(3.4));
  REQUIRE(b(4) == Approx(4.2));
  REQUIRE(b(5) == Approx(5.0));
  
  rowvec c = linspace<rowvec>(1,5,6);
  
  REQUIRE(c(0) == Approx(1.0));
  REQUIRE(c(1) == Approx(1.8));
  REQUIRE(c(2) == Approx(2.6));
  REQUIRE(c(3) == Approx(3.4));
  REQUIRE(c(4) == Approx(4.2));
  REQUIRE(c(5) == Approx(5.0));
  
  mat X = linspace<mat>(1,5,6);
  
  REQUIRE(X.n_rows == 6);
  REQUIRE(X.n_cols == 1);
  }



