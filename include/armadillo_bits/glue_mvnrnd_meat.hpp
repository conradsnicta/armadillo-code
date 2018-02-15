// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
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


//! \addtogroup glue_mvnrnd
//! @{


template<typename T1, typename T2>
inline
void
glue_mvnrnd::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_mvnrnd>& expr)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_mvnrnd::apply_direct_mode1(out, expr.A, expr.B, expr.aux_uword);
  
  if(status == false)
    {
    arma_stop_runtime_error("mvnrnd(): given covariance matrix is not symmetric positive semi-definite");
    }
  }



template<typename T1, typename T2>
inline
bool
glue_mvnrnd::apply_direct_mode1(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& M, const Base<typename T1::elem_type,T2>& C, const uword N)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UM(M.get_ref());
  const quasi_unwrap<T2> UC(C.get_ref());
  
  arma_debug_check( (UM.M.is_colvec() == false),  "mvnrnd(): given mean must be a column vector"                                   );
  arma_debug_check( (UC.M.is_square() == false),  "mvnrnd(): given covariance matrix must be square sized"                         );
  arma_debug_check( (UM.M.n_rows != UC.M.n_rows), "mvnrnd(): number of rows in given mean vector and covariance matrix must match" );
  
  bool status = false;
  
  if(UM.is_alias(out) || UC.is_alias(out))
    {
    Mat<eT> tmp;
    
    status = glue_mvnrnd::apply_noalias_mode1(tmp, UM.M, UC.M, N);
    
    out.steal_mem(tmp);
    }
  else
    {
    status = glue_mvnrnd::apply_noalias_mode1(out, UM.M, UC.M, N);
    }
  
  if(status == false)  { out.soft_reset(); }
  
  return status;
  }



template<typename T1, typename T2>
inline
bool
glue_mvnrnd::apply_direct_mode2(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& M, const Base<typename T1::elem_type,T2>& D, const uword N)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UM(M.get_ref());
  const quasi_unwrap<T2> UD(D.get_ref());
  
  arma_debug_check( (UM.M.is_colvec() == false),  "mvnrnd(): given mean must be a column vector"                                 );
  arma_debug_check( (UD.M.is_square() == false),  "mvnrnd(): given cholesky matrix must be square sized"                         );
  arma_debug_check( (UM.M.n_rows != UD.M.n_rows), "mvnrnd(): number of rows in given mean vector and cholesky matrix must match" );
  
  bool status = false;
  
  if(UM.is_alias(out) || UD.is_alias(out))
    {
    Mat<eT> tmp;
    
    status = glue_mvnrnd::apply_noalias_mode2(tmp, UM.M, UD.M, N);
    
    out.steal_mem(tmp);
    }
  else
    {
    status = glue_mvnrnd::apply_noalias_mode2(out, UM.M, UD.M, N);
    }
  
  if(status == false)  { out.soft_reset(); }
  
  return status;
  }



template<typename eT>
inline
bool
glue_mvnrnd::apply_noalias_mode1(Mat<eT>& out, const Mat<eT>& M, const Mat<eT>& C, const uword N)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> D;
  
  const bool chol_status = op_chol::apply_direct(D, C, 1);
  
  if(chol_status == false)
    {
    // C is not symmetric positive definite, so try approximation based on diagonalisation
    
    Col<eT> eigval;  // NOTE: eT is constrained to be real (ie. float or double) in fn_mvnrnd.hpp
    Mat<eT> eigvec;
    
    const bool eig_status = auxlib::eig_sym_dc(eigval, eigvec, C);
    
    if(eig_status == false)  { return false; }
    
          eT*   eigval_mem    = eigval.memptr();
    const uword eigval_n_elem = eigval.n_elem;
    
    // since we're doing an approximation, tolerate tiny negative eigenvalues
    
    const eT tol = eT(-100) * Datum<eT>::eps * norm(C, "fro");
    
    if(arma_isfinite(tol) == false)  { return false; }
    
    // cout << "*** eps: " << Datum<eT>::eps << endl;
    // cout << "*** tol: " << tol << endl;
    // eigval.print("*** eigval:");
    
    for(uword i=0; i<eigval_n_elem; ++i)
      {
      const eT val = eigval_mem[i];
      
      if( (val < tol) || (arma_isfinite(val) == false) )  { return false; }
      }
    
    for(uword i=0; i<eigval_n_elem; ++i)  { if(eigval_mem[i] < eT(0))  { eigval_mem[i] = eT(0); } }
    
    D = eigvec * diagmat(sqrt(eigval));
    
    // D.print("***D:");
    }
  
  return glue_mvnrnd::apply_noalias_mode2(out, M, D, N);
  }



template<typename eT>
inline
bool
glue_mvnrnd::apply_noalias_mode2(Mat<eT>& out, const Mat<eT>& M, const Mat<eT>& D, const uword N)
  {
  arma_extra_debug_sigprint();
  
  out = D * randn< Mat<eT> >(M.n_rows, N);
  
  if(N == 1)
    {
    out += M;
    }
  else
  if(N > 1)
    {
    out.each_col() += M;
    }
  
  return true;
  }



//! @}
