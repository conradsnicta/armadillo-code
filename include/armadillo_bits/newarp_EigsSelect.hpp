// Copyright (C) 2016 National ICT Australia (NICTA)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
//
// Written by Yixuan Qiu


namespace newarp
{


//! The enumeration of selection rules of desired eigenvalues.
struct EigsSelect
  {
  enum SELECT_EIGENVALUE
    {
    LARGEST_MAGN = 0,  //!< Select eigenvalues with largest magnitude.
                       //!< Magnitude means the absolute value for real numbers and norm for complex numbers.
                       //!< Applies to both symmetric and general eigen solvers.

    LARGEST_REAL,      //!< Select eigenvalues with largest real part. Only for general eigen solvers.

    LARGEST_IMAG,      //!< Select eigenvalues with largest imaginary part (in magnitude). Only for general eigen solvers.

    LARGEST_ALGE,      //!< Select eigenvalues with largest algebraic value, considering any negative sign. Only for symmetric eigen solvers.

    SMALLEST_MAGN,     //!< Select eigenvalues with smallest magnitude. Applies to both symmetric and general eigen solvers.

    SMALLEST_REAL,     //!< Select eigenvalues with smallest real part. Only for general eigen solvers.

    SMALLEST_IMAG,     //!< Select eigenvalues with smallest imaginary part (in magnitude). Only for general eigen solvers.

    SMALLEST_ALGE,     //!< Select eigenvalues with smallest algebraic value. Only for symmetric eigen solvers.

    BOTH_ENDS          //!< Select eigenvalues half from each end of the spectrum.
                       //!< When `nev` is odd, compute more from the high end. Only for symmetric eigen solvers.
    };
  };


}  // namespace newarp
