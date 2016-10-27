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


TEST_CASE("fn_find_nonfinite_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745   0.051408;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153   0.035437;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317  -0.454499;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040   0.373833;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768   0.258704;\
    ";
  
  mat B = A;
  
  B( 6) =  datum::nan;
  B( 8) =  datum::inf;
  B(10) = -datum::inf;
  
  uvec indices1 = find_nonfinite(A);
  uvec indices2 = find_nonfinite(B);
  
  REQUIRE( indices1.n_elem == 0 );
  REQUIRE( indices2.n_elem == 3 );
  
  REQUIRE( accu(indices2 - uvec({6,8,10})) == 0 );
  
  // REQUIRE_THROWS(  );
  }
