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


TEST_CASE("fn_cond_1")
  {
  mat A =
    {
    { -0.78838,  0.69298,  0.41084,  0.90142 },
    {  0.49345, -0.12020,  0.78987,  0.53124 },
    {  0.73573,  0.52104, -0.22263,  0.40163 }
    };
  
  REQUIRE( cond(A)                      == Approx(1.7455) );
  REQUIRE( cond(A(span::all,span::all)) == Approx(1.7455) );
  }



TEST_CASE("fn_cond_2")
  {
  mat A = zeros<mat>(5,6);
  
  REQUIRE( is_finite(cond(A)) == false );
  }
