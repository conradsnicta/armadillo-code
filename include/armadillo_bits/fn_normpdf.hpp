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


//! \addtogroup fn_normpdf
//! @{



template<typename eT>
arma_inline
typename enable_if2< (is_real<eT>::value), eT >::result
normpdf(const eT x)
  {
  const eT out = std::exp(-0.5 * (x*x)) / Datum<eT>::sqrt2pi;
  
  return out;
  }



template<typename eT>
inline
typename enable_if2< (is_real<eT>::value), eT >::result
normpdf(const eT x, const eT mu, const eT sigma)
  {
  const eT tmp = (x - mu) / sigma;
  
  const eT out = std::exp(-0.5 * (tmp*tmp)) / ( sigma * Datum<eT>::sqrt2pi );
  
  return out;
  }



template<typename T1>
inline
typename enable_if2< (is_real<typename T1::elem_type>::value), Mat<typename T1::elem_type> >::result
normpdf(const Base<typename T1::elem_type, T1>& X_expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> tmp1(X_expr.get_ref());
  
  const Mat<eT>& X = tmp1.M;
  
  Mat<eT> out( size(X) );
  
  eT* out_mem = out.memptr();
  
  const uword N     = X.n_elem;
  const eT*   X_mem = X.memptr();
  
  const bool use_mp = arma_config::cxx11 && arma_config::openmp && mp_gate<eT,true>::eval(N);
  
  if(use_mp)
    {
    #if defined(ARMA_USE_OPENMP)
      {
      const int n_threads = mp_thread_limit::get();
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword i=0; i<N; ++i)
        {
        const eT tmp = X_mem[i];
        
        out_mem[i] = std::exp(-0.5 * (tmp*tmp)) / Datum<eT>::sqrt2pi;
        }
      }
    #endif
    }
  else
    {
    for(uword i=0; i<N; ++i)
      {
      const eT tmp = X_mem[i];
      
      out_mem[i] = std::exp(-0.5 * (tmp*tmp)) / Datum<eT>::sqrt2pi;
      }
    }
  
  return out;
  }



template<typename T1, typename T2, typename T3>
inline
typename enable_if2< (is_real<typename T1::elem_type>::value), Mat<typename T1::elem_type> >::result
normpdf(const Base<typename T1::elem_type, T1>& X_expr, const Base<typename T1::elem_type, T2>& M_expr, const Base<typename T1::elem_type, T3>& S_expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> tmp1(X_expr.get_ref());
  const quasi_unwrap<T2> tmp2(M_expr.get_ref());
  const quasi_unwrap<T3> tmp3(S_expr.get_ref());
  
  const Mat<eT>& X = tmp1.M;
  const Mat<eT>& M = tmp2.M;
  const Mat<eT>& S = tmp3.M;
  
  arma_debug_check( ((size(X) != size(M)) || (size(M) != size(S))), "normpdf(): size mismatch" );

  Mat<eT> out( size(X) );
  
  eT* out_mem = out.memptr();
  
  const uword N     = X.n_elem;
  const eT*   X_mem = X.memptr();
  const eT*   M_mem = M.memptr();
  const eT*   S_mem = S.memptr();
  
  const bool use_mp = arma_config::cxx11 && arma_config::openmp && mp_gate<eT,true>::eval(N);
  
  if(use_mp)
    {
    #if defined(ARMA_USE_OPENMP)
      {
      const int n_threads = mp_thread_limit::get();
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword i=0; i<N; ++i)
        {
        const eT sigma = S_mem[i];
        
        const eT tmp = (X_mem[i] - M_mem[i]) / sigma;
        
        out_mem[i] = std::exp(-0.5 * (tmp*tmp)) / (sigma * Datum<eT>::sqrt2pi);
        }
      }
    #endif
    }
  else
    {
    for(uword i=0; i<N; ++i)
      {
      const eT sigma = S_mem[i];
      
      const eT tmp = (X_mem[i] - M_mem[i]) / sigma;
      
      out_mem[i] = std::exp(-0.5 * (tmp*tmp)) / (sigma * Datum<eT>::sqrt2pi);
      }
    }
  
  return out;
  }



//! @}
