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


TEST_CASE("fn_all_1")
  {
  vec a(5, fill::zeros);
  vec b(5, fill::zeros);  b(0) = 1.0;
  vec c(5, fill::ones );
  
  REQUIRE( all(a) == false);
  REQUIRE( all(b) == false);
  REQUIRE( all(c) == true );
  
  REQUIRE( all(a(span::all)) == false);
  REQUIRE( all(b(span::all)) == false);
  REQUIRE( all(c(span::all)) == true );
  
  REQUIRE( all(  c -  c) == false);
  REQUIRE( all(2*c -2*c) == false);
  
  REQUIRE( all(c < 0.5) == false);
  REQUIRE( all(c > 0.5) == true );
  }



TEST_CASE("fn_all_2")
  {
  mat A(5, 6, fill::zeros);
  mat B(5, 6, fill::zeros);  B(0,0) = 1.0;
  mat C(5, 6, fill::ones );
  
  REQUIRE( all(vectorise(A)) == false);
  REQUIRE( all(vectorise(B)) == false);
  REQUIRE( all(vectorise(C)) == true );
  
  REQUIRE( all(vectorise(A(span::all,span::all))) == false);
  REQUIRE( all(vectorise(B(span::all,span::all))) == false);
  REQUIRE( all(vectorise(C(span::all,span::all))) == true );


  REQUIRE( all(vectorise(  C -  C)) == false);
  REQUIRE( all(vectorise(2*C -2*C)) == false);
  
  REQUIRE( all(vectorise(C) < 0.5) == false);
  REQUIRE( all(vectorise(C) > 0.5) == true );
  }



TEST_CASE("fn_all_3")
  {
  mat A(5, 6, fill::zeros);
  mat B(5, 6, fill::zeros);  B(0,0) = 1.0;
  mat C(5, 6, fill::ones );
  mat D(5, 6, fill::ones );  D(0,0) = 0.0;
  
  REQUIRE( accu( all(A)   == urowvec({0, 0, 0, 0, 0, 0}) ) == 6 );
  REQUIRE( accu( all(A,0) == urowvec({0, 0, 0, 0, 0, 0}) ) == 6 );
  REQUIRE( accu( all(A,1) == uvec   ({0, 0, 0, 0, 0}   ) ) == 5 );

  REQUIRE( accu( all(B)   == urowvec({0, 0, 0, 0, 0, 0}) ) == 6 );
  REQUIRE( accu( all(B,0) == urowvec({0, 0, 0, 0, 0, 0}) ) == 6 );
  REQUIRE( accu( all(B,1) == uvec   ({0, 0, 0, 0, 0}   ) ) == 5 );
  
  REQUIRE( accu( all(C)   == urowvec({1, 1, 1, 1, 1, 1}) ) == 6 );
  REQUIRE( accu( all(C,0) == urowvec({1, 1, 1, 1, 1, 1}) ) == 6 );
  REQUIRE( accu( all(C,1) == uvec   ({1, 1, 1, 1, 1}   ) ) == 5 );
  
  REQUIRE( accu( all(D)   == urowvec({0, 1, 1, 1, 1, 1}) ) == 6 );
  REQUIRE( accu( all(D,0) == urowvec({0, 1, 1, 1, 1, 1}) ) == 6 );
  REQUIRE( accu( all(D,1) == uvec   ({0, 1, 1, 1, 1}   ) ) == 5 );
  }


