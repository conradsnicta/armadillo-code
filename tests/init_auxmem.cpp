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


TEST_CASE("init_auxmem_1")
  {
  double data[] = { 1, 2, 3, 4, 5, 6 };
  
  mat A(data, 2, 3);
  mat B(data, 2, 3, false);
  mat C(data, 2, 3, false, true);
  
  REQUIRE( A(0,0) == double(1) );
  REQUIRE( A(1,0) == double(2) );
  
  REQUIRE( A(0,1) == double(3) );
  REQUIRE( A(1,1) == double(4) );
  
  REQUIRE( A(0,2) == double(5) );
  REQUIRE( A(1,2) == double(6) );
  
  A(0,0) = 123.0;  REQUIRE( data[0] == 1 );
  
  B(0,0) = 123.0;  REQUIRE( data[0] == 123.0 );
  
  REQUIRE_THROWS( C.set_size(5,6) );
  }



