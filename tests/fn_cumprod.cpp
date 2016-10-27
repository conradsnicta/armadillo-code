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


TEST_CASE("fn_cumprod_1")
  {
  colvec a = linspace<colvec>(1,5,6);
  rowvec b = linspace<rowvec>(1,5,6);
  
  colvec c = { 1.0000, 1.8000, 4.6800, 15.9120, 66.8304, 334.1520 };
  
  REQUIRE( accu(abs(cumprod(a) - c    )) == Approx(0.0) );
  REQUIRE( accu(abs(cumprod(b) - c.t())) == Approx(0.0) );
  
  REQUIRE_THROWS( b = cumprod(a) );
  }



TEST_CASE("fn_cumprod_2")
  {
  mat A =
    {
    { -0.78838,  0.69298,  0.41084,  0.90142 },
    {  0.49345, -0.12020,  0.78987,  0.53124 },
    {  0.73573,  0.52104, -0.22263,  0.40163 }
    };
  
  mat B =
    {
    { -0.788380,  0.692980,  0.410840,  0.901420 },
    { -0.389026, -0.083296,  0.324510,  0.478870 },
    { -0.286218, -0.043401, -0.072246,  0.192329 }

    };
  
  mat C =
    {
    { -0.788380, -0.546332, -0.224455, -0.202328 },
    {  0.493450, -0.059313, -0.046849, -0.024888 },
    {  0.735730,  0.383345, -0.085344, -0.034277 }
    };
  
  REQUIRE( accu(abs(cumprod(A)   - B)) == Approx(0.0) );
  REQUIRE( accu(abs(cumprod(A,0) - B)) == Approx(0.0) );
  REQUIRE( accu(abs(cumprod(A,1) - C)) == Approx(0.0) );
  }



