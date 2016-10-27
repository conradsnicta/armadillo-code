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


TEST_CASE("fn_clamp_1")
  {
  mat A = randu<mat>(5,6);
  
  mat B = clamp(A, 0.2, 0.8); 
  REQUIRE( B.min() == Approx(0.2) );
  REQUIRE( B.max() == Approx(0.8) );
  
  mat C = clamp(A, A.min(), 0.8); 
  REQUIRE( C.min() == A.min()     );
  REQUIRE( C.max() == Approx(0.8) );
  
  mat D = clamp(A, 0.2, A.max());   
  REQUIRE( D.min() == Approx(0.2) );
  REQUIRE( D.max() == A.max()     );
  
  REQUIRE_THROWS( clamp(A, A.max(), A.min() ) );
  }
