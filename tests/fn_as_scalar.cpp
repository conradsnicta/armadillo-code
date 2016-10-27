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


TEST_CASE("fn_as_scalar_1")
  {
  mat A(1,1); A.fill(2.0);
  mat B(2,2); B.fill(2.0);
  
  REQUIRE( as_scalar(A) == Approx(2.0) );
  
  REQUIRE( as_scalar(2+A) == Approx(4.0) );
  
  REQUIRE( as_scalar(B(span(0,0), span(0,0))) == Approx(2.0) );
  
  REQUIRE_THROWS( as_scalar(B) );
  }



TEST_CASE("fn_as_scalar_2")
  {
  rowvec r = linspace<rowvec>(1,5,6);
  colvec q = linspace<colvec>(1,5,6);
  mat    X = 0.5*toeplitz(q);
  
  REQUIRE( as_scalar(r*q) == Approx(65.2) );
  
  REQUIRE( as_scalar(r*X*q) == Approx(380.848) );
  
  REQUIRE( as_scalar(r*diagmat(X)*q) == Approx(32.6) );
  REQUIRE( as_scalar(r*inv(diagmat(X))*q) == Approx(130.4) );
  }



TEST_CASE("fn_as_scalar_3")
  {
  cube A(1,1,1); A.fill(2.0);
  cube B(2,2,2); B.fill(2.0);
  
  REQUIRE( as_scalar(A) == Approx(2.0) );
  
  REQUIRE( as_scalar(2+A) == Approx(4.0) );
  
  REQUIRE( as_scalar(B(span(0,0), span(0,0), span(0,0))) == Approx(2.0) );
  
  REQUIRE_THROWS( as_scalar(B) );
  }
