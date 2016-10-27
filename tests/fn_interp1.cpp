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


TEST_CASE("fn_interp1_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745   0.051408;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153   0.035437;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317  -0.454499;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040   0.373833;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768   0.258704;\
    ";
  
  vec x = vectorise(A(0,0,size(5,3)));
  vec y = vectorise(A(0,3,size(5,3)));
  
  vec xi_a = linspace<vec>( x.min(), x.max(), 10 );
  vec xi_b = flipud( linspace<vec>( x.min(), x.max(), 11 ) );
  
  vec yi_a;
  vec yi_b;
  
  interp1(x, y, xi_a, yi_a);
  interp1(x, y, xi_b, yi_b);
  
  vec yi_a_gt = { 0.419733, 0.241248, 0.149666, 0.058084, 0.057588, 0.152062, -0.284524, -0.307613, -0.336627, 0.373833 };
  vec yi_b_gt = { 0.373833, -0.300357, -0.353940, -0.201854, -0.449865, 0.063571, 0.045817, 0.085559, 0.167982, 0.250406, 0.419733 };
  
  REQUIRE( accu(abs( yi_a - yi_a_gt )) == Approx(0.0) );
  REQUIRE( accu(abs( yi_b - yi_b_gt )) == Approx(0.0) );
  
  // REQUIRE_THROWS(  );
  }
