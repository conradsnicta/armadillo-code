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


TEST_CASE("init_misc_1")
  {
  const uword n_rows = 5;
  const uword n_cols = 6;
  
  mat A(5,6, fill::zeros);
  
  for(uword row=0; row<n_rows; ++row)
  for(uword col=0; col<n_cols; ++col)
    {
    REQUIRE(A(row,col) == 0);
    }
  
  mat B(5,6, fill::ones);
  
  for(uword row=0; row<n_rows; ++row)
  for(uword col=0; col<n_cols; ++col)
    {
    REQUIRE(B(row,col) == 1);
    }
  
  mat C(2,3);
  
  C(0,0) = 0;
  C(1,0) = 1;
  
  C(0,1) = 2;
  C(1,1) = 3;
  
  C(0,2) = 4;
  C(1,2) = 5;
  
  REQUIRE(C(0,0) == 0);
  REQUIRE(C(1,0) == 1);
  
  REQUIRE(C(0,1) == 2);
  REQUIRE(C(1,1) == 3);
  
  REQUIRE(C(0,2) == 4);
  REQUIRE(C(1,2) == 5);
  
  mat D = { {0, 2, 4},
            {1, 3, 5} };
  
  REQUIRE(D(0,0) == 0);
  REQUIRE(D(1,0) == 1);
  
  REQUIRE(D(0,1) == 2);
  REQUIRE(D(1,1) == 3);
  
  REQUIRE(D(0,2) == 4);
  REQUIRE(D(1,2) == 5);
  
  mat E;
  E << 0 << 2 << 4 << endr
    << 1 << 3 << 5 << endr;
  
  REQUIRE(E(0,0) == 0);
  REQUIRE(E(1,0) == 1);
  
  REQUIRE(E(0,1) == 2);
  REQUIRE(E(1,1) == 3);
  
  REQUIRE(E(0,2) == 4);
  REQUIRE(E(1,2) == 5);
  }



TEST_CASE("init_misc_2")
  {
  mat A =
    {
    { -0.78838,  0.69298,  0.41084,  0.90142 },
    {  0.49345, -0.12020,  0.78987,  0.53124 },
    {  0.73573,  0.52104, -0.22263,  0.40163 }
    };
  
  mat B =
    "\
    -0.78838,  0.69298,  0.41084,  0.90142;\
     0.49345, -0.12020,  0.78987,  0.53124;\
     0.73573,  0.52104, -0.22263,  0.40163;\
    ";
  
  mat C =
    "\
    -0.78838  0.69298  0.41084  0.90142;\
     0.49345 -0.12020  0.78987  0.53124;\
     0.73573  0.52104 -0.22263  0.40163;\
    ";
  
  mat D;
  D << -0.78838 <<  0.69298 <<  0.41084 << 0.90142 << endr
    <<  0.49345 << -0.12020 <<  0.78987 << 0.53124 << endr
    <<  0.73573 <<  0.52104 << -0.22263 << 0.40163 << endr;
  
  REQUIRE( A.n_rows == 3 );
  REQUIRE( A.n_cols == 4 );
  
  REQUIRE( A(0,0) == Approx(-0.78838) );
  REQUIRE( A(1,0) == Approx( 0.49345) );
  REQUIRE( A(2,0) == Approx( 0.73573) );
  REQUIRE( A(0,1) == Approx( 0.69298) );
  REQUIRE( A(1,1) == Approx(-0.12020) );
  REQUIRE( A(2,1) == Approx( 0.52104) );
  REQUIRE( A(0,2) == Approx( 0.41084) );
  REQUIRE( A(1,2) == Approx( 0.78987) );
  REQUIRE( A(2,2) == Approx(-0.22263) );
  REQUIRE( A(0,3) == Approx( 0.90142) );
  REQUIRE( A(1,3) == Approx( 0.53124) );
  REQUIRE( A(2,3) == Approx( 0.40163) );
  
  REQUIRE( accu(abs(A-B)) == Approx(0.0) );
  REQUIRE( accu(abs(A-C)) == Approx(0.0) );
  REQUIRE( accu(abs(A-D)) == Approx(0.0) );
  }



TEST_CASE("init_misc_3")
  {
  const uword n_rows = 5;
  const uword n_cols = 6;
  
  cx_mat A(5,6, fill::zeros);
  
  for(uword row=0; row<n_rows; ++row)
  for(uword col=0; col<n_cols; ++col)
    {
    REQUIRE(A(row,col) == cx_double(0,0));
    }
  
  cx_mat B(5,6, fill::ones);
  
  for(uword row=0; row<5; ++row)
  for(uword col=0; col<6; ++col)
    {
    REQUIRE(B(row,col) == cx_double(1,0));
    }
  
  cx_mat C(2,3);
  
  C(0,0) = cx_double(0.0,  1.0);
  C(1,0) = cx_double(1.0,  2.0);
  
  C(0,1) = cx_double(3.0,  4.0);
  C(1,1) = cx_double(5.0,  6.0);
  
  C(0,2) = cx_double(7.0,  8.0);
  C(1,2) = cx_double(9.0, 10.0);
  
  REQUIRE(C(0,0) == cx_double(0.0,  1.0));
  REQUIRE(C(1,0) == cx_double(1.0,  2.0));
  
  REQUIRE(C(0,1) == cx_double(3.0,  4.0));
  REQUIRE(C(1,1) == cx_double(5.0,  6.0));
  
  REQUIRE(C(0,2) == cx_double(7.0,  8.0));
  REQUIRE(C(1,2) == cx_double(9.0, 10.0));
  }



