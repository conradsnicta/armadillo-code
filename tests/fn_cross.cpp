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


TEST_CASE("fn_cross_1")
  {
  vec a = { 0.1,  2.3,  4.5 };
  vec b = { 6.7,  8.9, 10.0 };
  
  vec c = {-17.050, 29.150, -14.520 };
  
  REQUIRE( accu(abs(cross(a,b) - c)) == Approx(0.0) );
  
  vec x;
  
  REQUIRE_THROWS( x = cross(randu<vec>(4), randu<vec>(4)) );
  }



