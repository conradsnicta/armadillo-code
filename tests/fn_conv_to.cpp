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


TEST_CASE("fn_conv_to_1")
  {
  typedef std::vector<double> stdvec;

  stdvec x(3);
  x[0] = 10.0; x[1] = 20.0;  x[2] = 30.0;

  colvec y = conv_to< colvec >::from(x);
  stdvec z = conv_to< stdvec >::from(y);
  
  REQUIRE( z[0] == Approx(10.0) );
  REQUIRE( z[1] == Approx(20.0) );
  REQUIRE( z[2] == Approx(30.0) );
  }



TEST_CASE("fn_conv_to2")
  {
  mat A(5,6); A.fill(0.1);
  
  umat uA = conv_to<umat>::from(A);
  imat iA = conv_to<imat>::from(A);
  
  REQUIRE( (uA.n_rows - A.n_rows) == 0 );
  REQUIRE( (iA.n_rows - A.n_rows) == 0 );
  
  REQUIRE( (uA.n_cols - A.n_cols) == 0 );
  REQUIRE( (iA.n_cols - A.n_cols) == 0 );
  
  REQUIRE( any(vectorise(uA)) == false);
  REQUIRE( any(vectorise(iA)) == false);
  }
  

TEST_CASE("fn_conv_to3")
  {
  mat A(5,6); A.fill(1.0);
  
  umat uA = conv_to<umat>::from(A);
  imat iA = conv_to<imat>::from(A);
  
  REQUIRE( all(vectorise(uA)) == true);
  REQUIRE( all(vectorise(iA)) == true);
  }


TEST_CASE("fn_conv_to4")
  {
  mat A =   linspace<rowvec>(1,5,6);
  mat B = 2*linspace<colvec>(1,5,6);
  mat C = randu<mat>(5,6);
  
  REQUIRE( as_scalar( conv_to<rowvec>::from(A) * conv_to<colvec>::from(B) ) == Approx(130.40) );
  
  REQUIRE( conv_to<double>::from(A * B) == Approx(130.40) );
  
  REQUIRE_THROWS( conv_to<colvec>::from(C) );
  }


