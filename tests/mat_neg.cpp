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


TEST_CASE("mat_neg_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745   0.051408;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153   0.035437;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317  -0.454499;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040   0.373833;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768   0.258704;\
    ";
  
  mat B = -A;
  
  REQUIRE( B(0,0) == Approx(-0.061198) );
  REQUIRE( B(1,0) == Approx(-0.437242) );
  REQUIRE( B(2,0) == Approx(+0.492474) );
  REQUIRE( B(3,0) == Approx(-0.336352) );
  REQUIRE( B(4,0) == Approx(-0.239585) );
  
  REQUIRE( B(0,1) == Approx(-0.201990) );
  REQUIRE( B(1,1) == Approx(-0.058956) );
  REQUIRE( B(2,1) == Approx(+0.031309) );
  REQUIRE( B(3,1) == Approx(-0.411541) );
  REQUIRE( B(4,1) == Approx(+0.428913) );
  
  REQUIRE( B(0,5) == Approx(-0.051408) );
  REQUIRE( B(1,5) == Approx(-0.035437) );
  REQUIRE( B(2,5) == Approx(+0.454499) );
  REQUIRE( B(3,5) == Approx(-0.373833) );
  REQUIRE( B(4,5) == Approx(-0.258704) );
  
  REQUIRE( accu(B + A) == Approx(0.0) );
  
  // REQUIRE_THROWS(  );
  }
