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


TEST_CASE("bounds_1")
  {
  const uword n_rows = 5;
  const uword n_cols = 6;
  
  mat A(n_rows, n_cols, fill::zeros);
  
  REQUIRE_NOTHROW( A(n_rows-1,n_cols-1) = 0 );
  
  // out of bounds access will throw unless ARMA_NO_DEBUG is defined
  REQUIRE_THROWS( A(n_rows,n_cols) = 0 );
  }



