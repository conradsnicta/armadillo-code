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


TEST_CASE("fn_eps_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768;\
    ";
  
  // NOTE: this may not be portable,
  //       due to minor differences in representations of floating numbers between architectures
  
  mat B = 
    "\
     6.93889390390723e-18   2.77555756156289e-17   3.46944695195361e-18   5.55111512312578e-17   2.77555756156289e-17;\
     5.55111512312578e-17   6.93889390390723e-18   2.77555756156289e-17   6.93889390390723e-18   5.55111512312578e-17;\
     5.55111512312578e-17   6.93889390390723e-18   5.55111512312578e-17   5.55111512312578e-17   1.38777878078145e-17;\
     5.55111512312578e-17   5.55111512312578e-17   5.55111512312578e-17   5.55111512312578e-17   2.77555756156289e-17;\
     2.77555756156289e-17   5.55111512312578e-17   5.55111512312578e-17   5.55111512312578e-17   5.55111512312578e-17;\
    ";
  
  REQUIRE( accu(abs(eps(A) - B)) == Approx(0.0) );
  }
