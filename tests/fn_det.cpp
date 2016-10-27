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


TEST_CASE("fn_det_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745   0.051408;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153   0.035437;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317  -0.454499;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040   0.373833;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768   0.258704;\
    ";
  
  REQUIRE( det(A(0,0,size(0,0))) == Approx(+1.0               ) );
  REQUIRE( det(A(0,0,size(1,1))) == Approx(+0.0611980000000000) );
  REQUIRE( det(A(0,0,size(2,2))) == Approx(-0.0847105222920000) );
  REQUIRE( det(A(0,0,size(3,3))) == Approx(-0.0117387923199772) );
  REQUIRE( det(A(0,0,size(4,4))) == Approx(+0.0126070917169865) );
  REQUIRE( det(A(0,0,size(5,5))) == Approx(+0.0100409091117668) );
  
  REQUIRE_THROWS( det(A) );
  }



TEST_CASE("fn_det_2")
  {
  mat A = toeplitz(linspace(1,5,6));
  
  REQUIRE( det(A) == Approx(-31.45728) );
  
  
  mat B(6, 6, fill::zeros); B.diag() = linspace(1,5,6);
  
  REQUIRE( det(B) == Approx(334.152) );
  
  REQUIRE( det(diagmat(B)) == Approx(334.152) );
  
  
  mat C(5,6, fill::randu);
  
  REQUIRE_THROWS( det(C) );
  
  REQUIRE_THROWS( det(diagmat(C)) );
  }


TEST_CASE("fn_det_")
  {
  mat A = toeplitz(linspace(1,5,6));
  
  double val;
  double sign;
  
  log_det(val, sign, A);
  
  REQUIRE( val  == Approx(3.44863) );
  REQUIRE( sign == Approx(-1.0)    );
  
  REQUIRE( (std::exp(val)*sign) == Approx( det(A) ) );
  
  mat B(5,6, fill::randu);
  
  REQUIRE_THROWS( log_det(val, sign, B) );
  }

