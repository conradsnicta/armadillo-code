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


TEST_CASE("fn_find_unique_1")
  {
  mat A = 
    {
    {  1,  3,  5,  6,  7 },
    {  2,  4,  5,  7,  8 },
    {  3,  5,  5,  6,  9 },
    };
  
  uvec indices = find_unique(A);
  
  uvec indices2 = { 0, 1, 2, 4, 5, 9, 10, 13, 14 };
  
  REQUIRE( indices.n_elem == indices2.n_elem );
  
  bool same = true;
  
  for(uword i=0; i < indices.n_elem; ++i)
    {
    if(indices(i) != indices2(i))  { same = false; break; }
    }
  
  REQUIRE( same == true );
  
  vec unique_elem = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  
  REQUIRE( accu(abs( A.elem(indices) - unique_elem )) == Approx(0.0) );
  
  // REQUIRE_THROWS(  );
  }



TEST_CASE("fn_find_unique_2")
  {
  cx_mat A = 
    {
    { cx_double(1,-1), cx_double(3, 2), cx_double(5, 2), cx_double(6, 1), cx_double(7,-1) },
    { cx_double(2, 1), cx_double(4, 4), cx_double(5, 2), cx_double(7,-1), cx_double(8, 1) },
    { cx_double(3, 2), cx_double(5, 1), cx_double(5, 3), cx_double(6, 1), cx_double(9,-9) }
    };
  
  uvec indices = find_unique(A);
  
  uvec indices2 = { 0, 1, 2, 4, 5, 6, 8, 9, 10, 13, 14 };
  
  REQUIRE( indices.n_elem == indices2.n_elem );
  
  bool same = true;
  
  for(uword i=0; i < indices.n_elem; ++i)
    {
    if(indices(i) != indices2(i))  { same = false; break; }
    }
  
  REQUIRE( same == true );
  
  cx_vec unique_elem =
    {
    cx_double(1,-1), 
    cx_double(2, 1), 
    cx_double(3, 2), 
    cx_double(4, 4),
    cx_double(5, 1),
    cx_double(5, 2),
    cx_double(5, 3),
    cx_double(6, 1),
    cx_double(7,-1),
    cx_double(8, 1),
    cx_double(9,-9)
    };
  
  REQUIRE( accu(abs( A.elem(indices) - unique_elem )) == Approx(0.0) );
  
  // REQUIRE_THROWS(  );
  }
