// Copyright 2011-2017 Ryan Curtin (http://www.ratml.org/)
// Copyright 2017 National ICT Australia (NICTA)
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

TEST_CASE("subview_min_test")
  {
  // We will assume subview.at() works and returns points within the bounds of
  // the matrix, so we just have to ensure the results are the same as
  // Mat.min()...
  for (size_t r = 50; r < 150; ++r)
    {
    mat x;
    x.randu(r, r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;
    uword x_subview_min3;

    const double mval = x.min(x_min);
    const double mval1 = x.submat(0, 0, r - 1, r - 1).min(x_subview_min1);
    const double mval2 = x.cols(0, r - 1).min(x_subview_min2);
    const double mval3 = x.rows(0, r - 1).min(x_subview_min3);

    REQUIRE( x_min == x_subview_min1 );
    REQUIRE( x_min == x_subview_min2 );
    REQUIRE( x_min == x_subview_min3 );

    REQUIRE( mval == Approx(mval1) );
    REQUIRE( mval == Approx(mval2) );
    REQUIRE( mval == Approx(mval3) );
    }
  }



TEST_CASE("subview_col_min_test")
  {
  for (size_t r = 10; r < 50; ++r)
    {
    vec x;
    x.randu(r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;

    const double mval = x.min(x_min);
    const double mval1 = x.submat(0, 0, r - 1, 0).min(x_subview_min1);
    const double mval2 = x.rows(0, r - 1).min(x_subview_min2);

    REQUIRE( x_min == x_subview_min1 );
    REQUIRE( x_min == x_subview_min2 );

    REQUIRE( mval == Approx(mval1) );
    REQUIRE( mval == Approx(mval2) );
    }
  }



TEST_CASE("subview_row_min_test")
  {
  for (size_t r = 10; r < 50; ++r)
    {
    rowvec x;
    x.randu(r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;

    const double mval = x.min(x_min);
    const double mval1 = x.submat(0, 0, 0, r - 1).min(x_subview_min1);
    const double mval2 = x.cols(0, r - 1).min(x_subview_min2);

    REQUIRE( x_min == x_subview_min1 );
    REQUIRE( x_min == x_subview_min2 );

    REQUIRE( mval == Approx(mval1) );
    REQUIRE( mval == Approx(mval2) );
    }
  }



TEST_CASE("subview_incomplete_min_test")
  {
  for (size_t r = 50; r < 150; ++r)
    {
    mat x;
    x.randu(r, r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;
    uword x_subview_min3;

    const double mval = x.min(x_min);
    const double mval1 = x.submat(1, 1, r - 2, r - 2).min(x_subview_min1);
    const double mval2 = x.cols(1, r - 2).min(x_subview_min2);
    const double mval3 = x.rows(1, r - 2).min(x_subview_min3);

    uword row, col;
    x.min(row, col);

    if (row != 0 && row != r - 1 && col != 0 && col != r - 1)
      {
      uword srow, scol;

      srow = x_subview_min1 % (r - 2);
      scol = x_subview_min1 / (r - 2);
      REQUIRE( x_min == (srow + 1) + r * (scol + 1) );
      REQUIRE( x_min == x_subview_min2 + r );

      srow = x_subview_min3 % (r - 2);
      scol = x_subview_min3 / (r - 2);
      REQUIRE( x_min == (srow + 1) + r * scol );

      REQUIRE( mval == Approx(mval1) );
      REQUIRE( mval == Approx(mval2) );
      REQUIRE( mval == Approx(mval3) );
      }
    }
  }



TEST_CASE("subview_incomplete_col_min_test")
  {
  for (size_t r = 10; r < 50; ++r)
    {
    vec x;
    x.randu(r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;

    const double mval = x.min(x_min);
    const double mval1 = x.submat(1, 0, r - 2, 0).min(x_subview_min1);
    const double mval2 = x.rows(1, r - 2).min(x_subview_min2);

    if (x_min != 0 && x_min != r - 1)
      {
      REQUIRE( x_min == x_subview_min1 + 1 );
      REQUIRE( x_min == x_subview_min2 + 1 );

      REQUIRE( mval == Approx(mval1) );
      REQUIRE( mval == Approx(mval2) );
      }
    }
  }



TEST_CASE("subview_incomplete_row_min_test")
  {
  for (size_t r = 10; r < 50; ++r)
    {
    rowvec x;
    x.randu(r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;

    const double mval = x.min(x_min);
    const double mval1 = x.submat(0, 1, 0, r - 2).min(x_subview_min1);
    const double mval2 = x.cols(1, r - 2).min(x_subview_min2);

    if (x_min != 0 && x_min != r - 1)
      {
      REQUIRE( x_min == x_subview_min1 + 1 );
      REQUIRE( x_min == x_subview_min2 + 1 );

      REQUIRE( mval == Approx(mval1) );
      REQUIRE( mval == Approx(mval2) );
      }
    }
  }



TEST_CASE("subview_cx_min_test")
  {
  // We will assume subview.at() works and returns points within the bounds of
  // the matrix, so we just have to ensure the results are the same as
  // Mat.min()...
  for (size_t r = 50; r < 150; ++r)
    {
    cx_mat x;
    x.randu(r, r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;
    uword x_subview_min3;

    const std::complex<double> mval = x.min(x_min);
    const std::complex<double> mval1 = x.submat(0, 0, r - 1, r - 1).min(x_subview_min1);
    const std::complex<double> mval2 = x.cols(0, r - 1).min(x_subview_min2);
    const std::complex<double> mval3 = x.rows(0, r - 1).min(x_subview_min3);

    REQUIRE( x_min == x_subview_min1 );
    REQUIRE( x_min == x_subview_min2 );
    REQUIRE( x_min == x_subview_min3 );

    REQUIRE( mval.real() == Approx(mval1.real()) );
    REQUIRE( mval.imag() == Approx(mval1.imag()) );
    REQUIRE( mval.real() == Approx(mval2.real()) );
    REQUIRE( mval.imag() == Approx(mval2.imag()) );
    REQUIRE( mval.real() == Approx(mval3.real()) );
    REQUIRE( mval.imag() == Approx(mval3.imag()) );
    }
  }



TEST_CASE("subview_cx_col_min_test")
  {
  for (size_t r = 10; r < 50; ++r)
    {
    cx_vec x;
    x.randu(r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;

    const std::complex<double> mval = x.min(x_min);
    const std::complex<double> mval1 = x.submat(0, 0, r - 1, 0).min(x_subview_min1);
    const std::complex<double> mval2 = x.rows(0, r - 1).min(x_subview_min2);

    REQUIRE( x_min == x_subview_min1 );
    REQUIRE( x_min == x_subview_min2 );

    REQUIRE( mval.real() == Approx(mval1.real()) );
    REQUIRE( mval.imag() == Approx(mval1.imag()) );
    REQUIRE( mval.real() == Approx(mval2.real()) );
    REQUIRE( mval.imag() == Approx(mval2.imag()) );
    }
  }



TEST_CASE("subview_cx_row_min_test")
  {
  for (size_t r = 10; r < 50; ++r)
    {
    cx_rowvec x;
    x.randu(r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;

    const std::complex<double> mval = x.min(x_min);
    const std::complex<double> mval1 = x.submat(0, 0, 0, r - 1).min(x_subview_min1);
    const std::complex<double> mval2 = x.cols(0, r - 1).min(x_subview_min2);

    REQUIRE( x_min == x_subview_min1 );
    REQUIRE( x_min == x_subview_min2 );

    REQUIRE( mval.real() == Approx(mval1.real()) );
    REQUIRE( mval.imag() == Approx(mval1.imag()) );
    REQUIRE( mval.real() == Approx(mval2.real()) );
    REQUIRE( mval.imag() == Approx(mval2.imag()) );
    }
  }



TEST_CASE("subview_cx_incomplete_min_test")
  {
  for (size_t r = 50; r < 150; ++r)
    {
    cx_mat x;
    x.randu(r, r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;
    uword x_subview_min3;

    const std::complex<double> mval = x.min(x_min);
    const std::complex<double> mval1 = x.submat(1, 1, r - 2, r - 2).min(x_subview_min1);
    const std::complex<double> mval2 = x.cols(1, r - 2).min(x_subview_min2);
    const std::complex<double> mval3 = x.rows(1, r - 2).min(x_subview_min3);

    uword row, col;
    x.min(row, col);

    if (row != 0 && row != r - 1 && col != 0 && col != r - 1)
      {
      uword srow, scol;

      srow = x_subview_min1 % (r - 2);
      scol = x_subview_min1 / (r - 2);
      REQUIRE( x_min == (srow + 1) + r * (scol + 1) );
      REQUIRE( x_min == x_subview_min2 + r );

      srow = x_subview_min3 % (r - 2);
      scol = x_subview_min3 / (r - 2);
      REQUIRE( x_min == (srow + 1) + r * scol );

      REQUIRE( mval.real() == Approx(mval1.real()) );
      REQUIRE( mval.imag() == Approx(mval1.imag()) );
      REQUIRE( mval.real() == Approx(mval2.real()) );
      REQUIRE( mval.imag() == Approx(mval2.imag()) );
      REQUIRE( mval.real() == Approx(mval3.real()) );
      REQUIRE( mval.imag() == Approx(mval3.imag()) );
      }
    }
  }



TEST_CASE("subview_cx_incomplete_col_min_test")
  {
  for (size_t r = 10; r < 50; ++r)
    {
    cx_vec x;
    x.randu(r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;

    const std::complex<double> mval = x.min(x_min);
    const std::complex<double> mval1 = x.submat(1, 0, r - 2, 0).min(x_subview_min1);
    const std::complex<double> mval2 = x.rows(1, r - 2).min(x_subview_min2);

    if (x_min != 0 && x_min != r - 1)
      {
      REQUIRE( x_min == x_subview_min1 + 1 );
      REQUIRE( x_min == x_subview_min2 + 1 );

      REQUIRE( mval.real() == Approx(mval1.real()) );
      REQUIRE( mval.imag() == Approx(mval1.imag()) );
      REQUIRE( mval.real() == Approx(mval2.real()) );
      REQUIRE( mval.imag() == Approx(mval2.imag()) );
      }
    }
  }



TEST_CASE("subview_cx_incomplete_row_min_test")
  {
  for (size_t r = 10; r < 50; ++r)
    {
    cx_rowvec x;
    x.randu(r);

    uword x_min;
    uword x_subview_min1;
    uword x_subview_min2;

    const std::complex<double> mval = x.min(x_min);
    const std::complex<double> mval1 = x.submat(0, 1, 0, r - 2).min(x_subview_min1);
    const std::complex<double> mval2 = x.cols(1, r - 2).min(x_subview_min2);

    if (x_min != 0 && x_min != r - 1)
      {
      REQUIRE( x_min == x_subview_min1 + 1 );
      REQUIRE( x_min == x_subview_min2 + 1 );

      REQUIRE( mval.real() == Approx(mval1.real()) );
      REQUIRE( mval.imag() == Approx(mval1.imag()) );
      REQUIRE( mval.real() == Approx(mval2.real()) );
      REQUIRE( mval.imag() == Approx(mval2.imag()) );
      }
    }
  }
