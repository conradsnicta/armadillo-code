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


TEST_CASE("gen_randu_1")
  {
  const uword n_rows = 100;
  const uword n_cols = 101;
  
  mat A(n_rows,n_cols, fill::randu);
  
  mat B(n_rows,n_cols); B.randu();
  
  mat C; C.randu(n_rows,n_cols);
  
  REQUIRE( (accu(A)/A.n_elem) == Approx(0.5).epsilon(0.01) );
  REQUIRE( (accu(B)/A.n_elem) == Approx(0.5).epsilon(0.01) );
  REQUIRE( (accu(C)/A.n_elem) == Approx(0.5).epsilon(0.01) );

  REQUIRE( (mean(vectorise(A))) == Approx(0.5).epsilon(0.01) );
  }
  


TEST_CASE("gen_randu_2")
  {
  mat A(50,60,fill::zeros);
  
  A(span(1,48),span(1,58)).randu();
  
  REQUIRE( accu(A.head_cols(1)) == Approx(0.0) );
  REQUIRE( accu(A.head_rows(1)) == Approx(0.0) );
  
  REQUIRE( accu(A.tail_cols(1)) == Approx(0.0) );
  REQUIRE( accu(A.tail_rows(1)) == Approx(0.0) );
  
  REQUIRE( mean(vectorise(A(span(1,48),span(1,58)))) == Approx(double(0.5)).epsilon(0.01) );
  }



