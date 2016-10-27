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


TEST_CASE("fn_conv_1")
  {
  vec a =   linspace<vec>(1,5,6);
  vec b = 2*linspace<vec>(1,6,7);
  
  vec c = conv(a,b);
  vec d =
    {
      2.00000000000000,
      7.26666666666667,
     17.13333333333333,
     32.93333333333334,
     56.00000000000000,
     87.66666666666667,
    117.66666666666666,
    134.00000000000003,
    137.73333333333335,
    127.53333333333336,
    102.06666666666668,
     60.00000000000000
     };
  
  REQUIRE( accu(abs(c - d)) == Approx(0.0) );
  }
