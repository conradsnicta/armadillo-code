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


TEST_CASE("fn_vectorise_1")
  {
  mat A = 
    "\
     0.061198   0.201990;\
     0.437242   0.058956;\
    -0.492474  -0.031309;\
     0.336352   0.411541;\
    ";
  
  vec a =
    {
     0.061198,
     0.437242,
    -0.492474,
     0.336352,
     0.201990,
     0.058956,
    -0.031309,
     0.411541,
    };
  
  rowvec b = { 0.061198, 0.201990, 0.437242, 0.058956, -0.492474, -0.031309, 0.336352, 0.411541 };
  
  REQUIRE( accu(abs(a - vectorise(A  ))) == Approx(0.0) );
  REQUIRE( accu(abs(a - vectorise(A,0))) == Approx(0.0) );
  REQUIRE( accu(abs(b - vectorise(A,1))) == Approx(0.0) );
  }
