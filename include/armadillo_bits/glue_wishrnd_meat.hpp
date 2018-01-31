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


//! \addtogroup glue_wishrnd
//! @{


template<typename T1, typename T2>
inline
void
glue_wishrnd::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_wishrnd>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword mode  = expr.aux_uword;
  const eT    df = expr.aux;
  
  if(mode == 1)
    {
    const quasi_unwrap<T1> UA(expr.A);
    
    if(UA.is_alias(out))
      {
      Mat<eT> tmp;
      
      glue_wishrnd::apply_noalias(tmp, UA.M, df);
      
      out.steal_mem(tmp);
      }
    else
      {
      glue_wishrnd::apply_noalias(out, UA.M, df);
      }
    }
  else
  if(mode == 2)
    {
    const quasi_unwrap<T1> UA(expr.A);
    const quasi_unwrap<T2> UB(expr.B);
    
    if(UA.is_alias(out) || UB.is_alias(out))
      {
      Mat<eT> tmp;
      
      glue_wishrnd::apply_noalias(tmp, UA.M, df, UB.M);
      
      out.steal_mem(tmp);
      }
    else
      {
      glue_wishrnd::apply_noalias(out, UA.M, df, UB.M);
      }
    }
  }


template<typename eT>
inline
void
glue_wishrnd::apply_noalias(Mat<eT>& out, const Mat<eT>& S, const eT df)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (S.is_square() == false), "wishrnd(): given matrix must be square sized" );
  
  Mat<eT> D;
  
  const bool status = op_chol::apply_direct(D, S, 0);
  
  // TODO: overload user-facing function to return bool on success
  // TODO: pattern after expmat() ?
  
  if(status == false)
    {
    out.reset();
    arma_debug_check(true, "wishrnd(): matrix is not positive definite");
    }
  
  glue_wishrnd::apply_noalias(out, S, df, D);
  }



template<typename eT>
inline
void
glue_wishrnd::apply_noalias(Mat<eT>& out, const Mat<eT>& S, const eT df, const Mat<eT>& D)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_CXX11)
    {
    arma_debug_check( (S.is_square() == false), "wishrnd(): given matrix must be square sized" );
    arma_debug_check( (size(S) != size(D)),     "wishrnd(): size mismatch"                     );
    
    if(S.is_empty())  { out.reset(); return; }
    
    const uword N = S.n_rows;
    
    const eT   df_floor  = std::floor(df);
    const eT   df_val    = (df < eT(N)) ? df_floor : df;
    const bool df_simple = (df == df_floor);
    
    if(df_simple)
      {
      arma_extra_debug_print("simple generator");
      
      const Mat<eT> tmp = (randn< Mat<eT> >(uword(df_val), N)) * D;
      
      out = tmp.t() * tmp;
      }
    else
      {
      arma_extra_debug_print("standard generator");
      
      op_chi2rnd_generator<eT> chi2rnd_gen;
      
      Mat<eT> A(N, N, fill::zeros);
      
      for(uword i=0; i<N; ++i)
        {
        A.at(i,i) = std::sqrt( chi2rnd_gen(df_val - eT(i)) );
        }
      
      for(uword j=  0; j<N; ++j)
      for(uword i=j+1; i<N; ++i)
        {
        A.at(i,j) = eT( arma_rng::randn<eT>() );
        }
      
      const Mat<eT> B = D * A;
      
      A.reset();
      
      out = B * B.t();
      }
    }
  #else
    {
    arma_ignore(S);
    arma_ignore(df);
    arma_ignore(D);
    arma_stop_logic_error("C++11 compiler required");
    out.reset();
    }
  #endif
  }



//! @}
