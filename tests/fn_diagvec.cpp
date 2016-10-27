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


TEST_CASE("fn_diagvec_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768;\
    ";
  
  vec A_main1 = diagvec(A);
  vec A_main2 = diagvec(A,0);
  
  vec A_p1 = diagvec(A, 1);
  vec A_m1 = diagvec(A,-1);
  
  vec a =
    {
     0.061198,
     0.058956,
     0.314156,
    -0.393139,
    -0.353768
    };
  
  vec b =
    {
     0.20199,
    -0.14936,
     0.41973,
    -0.13504
    };
  
  vec c =
    {
     0.437242,
    -0.031309,
     0.458476,
    -0.291020
    };
  
  REQUIRE( accu(abs(A_main1 - a)) == Approx(0.0) );
  REQUIRE( accu(abs(A_main2 - a)) == Approx(0.0) );
  
  REQUIRE( accu(abs(A_p1 - b)) == Approx(0.0) );
  REQUIRE( accu(abs(A_m1 - c)) == Approx(0.0) );
  }
