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


TEST_CASE("init_fill_1")
  {
  mat Z( 5,  6, fill::zeros);
  mat O( 5,  6, fill::ones);
  mat I( 5,  6, fill::eye);
  mat U(50, 60, fill::randu);
  mat N(50, 60, fill::randn);
  
  
  REQUIRE( accu(Z != 0) == 0   );
  REQUIRE( accu(O != 0) == 5*6 );
  REQUIRE( accu(I != 0) == 5   );
  
  REQUIRE(   mean(vectorise(U)) == Approx(0.500).epsilon(0.05) );
  REQUIRE( stddev(vectorise(U)) == Approx(0.288).epsilon(0.05) );
  
  REQUIRE(   mean(vectorise(N)) == Approx(0.0).epsilon(0.05) );
  REQUIRE( stddev(vectorise(N)) == Approx(1.0).epsilon(0.05) );
  
  mat X(5, 6, fill::none);   // only to test instantiation
  }



TEST_CASE("init_fill_2")
  {
  cube Z( 5,  6, 2, fill::zeros);
  cube O( 5,  6, 2, fill::ones);
  cube U(50, 60, 2, fill::randu);
  cube N(50, 60, 2, fill::randn);
  
  REQUIRE( accu(Z != 0) == 0     );
  REQUIRE( accu(O != 0) == 5*6*2 );
  
  REQUIRE(   mean(vectorise(U)) == Approx(0.500).epsilon(0.05) );
  REQUIRE( stddev(vectorise(U)) == Approx(0.288).epsilon(0.05) );
  
  REQUIRE(   mean(vectorise(N)) == Approx(0.0).epsilon(0.05) );
  REQUIRE( stddev(vectorise(N)) == Approx(1.0).epsilon(0.05) );
  
  cube X(5, 6, 2, fill::none);   // only to test instantiation
  
  cube I;  REQUIRE_THROWS( I = cube(5, 6, 2, fill::eye) );
  }
