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


TEST_CASE("fn_cov_1")
  {
  vec a =     linspace<vec>(1,5,6);
  vec b = 0.5*linspace<vec>(1,5,6);
  vec c = flipud(b);
  
  REQUIRE( as_scalar(cov(a,b) - (+1.12)) == Approx(0.0) );
  REQUIRE( as_scalar(cov(a,c) - (-1.12)) == Approx(0.0) );
  }



TEST_CASE("fn_cov_2")
  {
  mat A =
    {
    { -0.78838,  0.69298,  0.41084,  0.90142 },
    {  0.49345, -0.12020,  0.78987,  0.53124 },
    {  0.73573,  0.52104, -0.22263,  0.40163 }
    };
  
  mat B = 0.5 * A;

  mat C = fliplr(B);
  
  mat AA =
    "\
     0.670783  -0.191509  -0.120822  -0.211274;\
    -0.191509   0.183669  -0.141426   0.050641;\
    -0.120822  -0.141426   0.261684   0.051254;\
    -0.211274   0.050641   0.051254   0.067270;\
    ";
  
  mat AB =
    "\
     0.335392  -0.095755  -0.060411  -0.105637;\
    -0.095755   0.091834  -0.070713   0.025320;\
    -0.060411  -0.070713   0.130842   0.025627;\
    -0.105637   0.025320   0.025627   0.033635;\
    ";
    
  mat AC =
    "\
    -0.105637  -0.060411  -0.095755   0.335392;\
     0.025320  -0.070713   0.091834  -0.095755;\
     0.025627   0.130842  -0.070713  -0.060411;\
     0.033635   0.025627   0.025320  -0.105637;\
    ";
    
  REQUIRE( accu(abs(cov(A)   - AA)) == Approx(0.0) );
  REQUIRE( accu(abs(cov(A,B) - AB)) == Approx(0.0) );
  REQUIRE( accu(abs(cov(A,C) - AC)) == Approx(0.0) );
  }
