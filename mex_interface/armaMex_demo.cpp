// Copyright (C) 2014 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by George Yammine


// Demonstration of how to connect Armadillo with Matlab mex functions.
// Version 0.2


#include "armaMex.hpp"


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
  {
  // Check the number of input arguments.
  if (nrhs != 2)
    mexErrMsgTxt("Incorrect number of input arguments.");
  
  // Check type of input.
  if ( (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) || (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) )
    mexErrMsgTxt("Input must me of type double.");
  
  // Check if input is real.
  if ( (mxIsComplex(prhs[0])) || (mxIsComplex(prhs[1])) )
    mexErrMsgTxt("Input must be real.");
  
  // Create matrices X and Y from the first and second argument.
  mat X = armaGetPr(prhs[0]);
  mat Y = armaGetPr(prhs[1]);
  
  // Our calculations require that matrices must be of the same size 
  if ( size(X) != size(Y) )
    mexErrMsgTxt("Matrices should be of same size.");
  
  // Perform calculations
  mat A = X + Y;
  mat B = X % Y;  // % means element-wise multiplication in Armadillo
  
  // Create cube C with A and B as slices.
  cube C(A.n_rows, A.n_cols, 2);
  
  C.slice(0) = A;
  C.slice(1) = B;
  
  // Create the output argument plhs[0] to return cube C
  plhs[0] = armaCreateMxMatrix(C.n_rows, C.n_cols, C.n_slices);
  
  // Return the cube C as plhs[0] in Matlab/Octave
  armaSetCubePr(plhs[0], C);
  
  return;
  }
