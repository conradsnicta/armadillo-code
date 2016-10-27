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


TEST_CASE("fn_is_finite_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745   0.051408;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153   0.035437;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317  -0.454499;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040   0.373833;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768   0.258704;\
    ";
  
  mat B = A;  B(1,1) = datum::inf;
  
  mat C = A;  C(2,4) = datum::nan;
  
  REQUIRE( is_finite(A) == true  );
  REQUIRE( is_finite(B) == false );
  REQUIRE( is_finite(C) == false );
  
  REQUIRE( is_finite(A+A) == true  );
  REQUIRE( is_finite(B+B) == false );
  REQUIRE( is_finite(C+C) == false );
  
  REQUIRE( is_finite(2*A) == true  );
  REQUIRE( is_finite(2*B) == false );
  REQUIRE( is_finite(2*C) == false );
  
  // REQUIRE_THROWS(  );
  }
