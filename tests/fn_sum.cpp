// Copyright 2015 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2015 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


#include <armadillo>
#include "catch.hpp"

using namespace arma;


TEST_CASE("fn_sum_1")
  {
  vec a = linspace<vec>(1,5,5);
  vec b = linspace<vec>(1,5,6);
  
  REQUIRE(sum(a) == Approx(15.0));
  REQUIRE(sum(b) == Approx(18.0));
  }



TEST_CASE("sum2")
  {
  mat A =
    {
    { -0.78838,  0.69298,  0.41084,  0.90142 },
    {  0.49345, -0.12020,  0.78987,  0.53124 },
    {  0.73573,  0.52104, -0.22263,  0.40163 }
    };
  
  rowvec colsums = { 0.44080, 1.09382, 0.97808, 1.83429 };
  
  colvec rowsums =
    {
    1.21686,
    1.69436,
    1.43577
    };
  
  REQUIRE( accu(abs(colsums - sum(A  ))) == Approx(0.0) );
  REQUIRE( accu(abs(colsums - sum(A,0))) == Approx(0.0) );
  REQUIRE( accu(abs(rowsums - sum(A,1))) == Approx(0.0) );
  }


TEST_CASE("sum3")
  {
  mat AA =
    {
    { -0.78838,  0.69298,  0.41084,  0.90142 },
    {  0.49345, -0.12020,  0.78987,  0.53124 },
    {  0.73573,  0.52104, -0.22263,  0.40163 }
    };
  
  cx_mat A = cx_mat(AA, 0.5*AA);
  
  rowvec re_colsums = { 0.44080, 1.09382, 0.97808, 1.83429 };
  
  cx_rowvec cx_colsums = cx_rowvec(re_colsums, 0.5*re_colsums);
  
  colvec re_rowsums =
    {
    1.21686,
    1.69436,
    1.43577
    };
  
  cx_colvec cx_rowsums = cx_colvec(re_rowsums, 0.5*re_rowsums);
  
  REQUIRE( accu(abs(cx_colsums - sum(A  ))) == Approx(0.0) );
  REQUIRE( accu(abs(cx_colsums - sum(A,0))) == Approx(0.0) );
  REQUIRE( accu(abs(cx_rowsums - sum(A,1))) == Approx(0.0) );
  }


TEST_CASE("sum4")
  {
  mat X(100,101, fill::randu);
  
  REQUIRE( (sum(sum(X))/X.n_elem)                      == Approx(0.5).epsilon(0.01) );
  REQUIRE( (sum(sum(X(span::all,span::all)))/X.n_elem) == Approx(0.5).epsilon(0.01) );
  }
