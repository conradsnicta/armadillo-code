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


TEST_CASE("fn_dot_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768;\
    ";
  
  vec a = A.head_cols(1);
  vec b = A.tail_cols(1);
  
  rowvec c = A.head_rows(1);
  rowvec d = A.tail_rows(1);
  
  REQUIRE( dot(  a,  b) == Approx(-0.04208883710200) );
  REQUIRE( dot(2*a,2+b) == Approx( 2.24343432579600) );
  
  REQUIRE( dot(    c,  d) == Approx( 0.108601544706000) );
  REQUIRE( dot(0.5*c,2-d) == Approx(-0.392115772353000) );
  
  REQUIRE( dot(a,b) == Approx( dot(A.head_cols(1), A.tail_cols(1)) ) );
  REQUIRE( dot(c,d) == Approx( dot(A.head_rows(1), A.tail_rows(1)) ) );
  }



TEST_CASE("fn_dot_2")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768;\
    ";
  
  cx_vec a = cx_vec(A.col(0), A.col(1));
  cx_vec b = cx_vec(A.col(2), A.col(3));
  
  cx_rowvec c = cx_rowvec(A.row(0), A.row(1));
  cx_rowvec d = cx_rowvec(A.row(2), A.row(3));
  
  REQUIRE( abs( dot(a,b) - cx_double(-0.009544718641000, -0.110209641379000)) == Approx(0.0) );
  REQUIRE( abs( dot(c,d) - cx_double(-0.326993347830000, +0.061084261990000)) == Approx(0.0) );
  
  REQUIRE( abs(cdot(a,b) - cx_double(-0.314669805873000, -0.807333974477000)) == Approx(0.0) );
  REQUIRE( abs(cdot(c,d) - cx_double(-0.165527940664000, +0.586984291846000)) == Approx(0.0) );
  }


// TODO: norm_dot
