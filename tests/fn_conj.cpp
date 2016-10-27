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


TEST_CASE("fn_conj_1")
  {
  vec re =   linspace<vec>(1,5,6);
  vec im = 2*linspace<vec>(1,5,6);
  
  cx_vec a = cx_vec(re,im);
  cx_vec b = conj(a);
  
  REQUIRE( accu(abs(real(b) - ( re))) == Approx(0.0) );
  REQUIRE( accu(abs(imag(b) - (-im))) == Approx(0.0) );
  }



TEST_CASE("fn_conj2")
  {
  cx_mat A = randu<cx_mat>(5,6);
  
  cx_mat B = conj(A);
  
  REQUIRE( all(vectorise(real(B) ==  real(A))) == true );
  REQUIRE( all(vectorise(imag(B) == -imag(A))) == true );
  }
