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


//! \addtogroup fn_chi2rnd
//! @{



arma_warn_unused
inline
double
chi2rnd(const double df)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_CXX11)
    {
    op_chi2rnd_generator<double> generator;
    
    return generator(df);
    }
  #else
    {
    arma_stop_logic_error("chi2rnd(): C++11 compiler required");
    
    return double(0);
    }
  #endif
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_real<typename T1::elem_type>::value),
  const Op<T1, op_chi2rnd>
  >::result
chi2rnd(const T1& expr)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_chi2rnd>(expr);
  }



template<typename obj_type>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_Mat<obj_type>::value && is_real<typename obj_type::elem_type>::value),
  obj_type
  >::result
chi2rnd(const typename obj_type::elem_type df, const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_CXX11)
    {
    typedef typename obj_type::elem_type eT;
    
    return randg<obj_type>(n_rows, n_cols, distr_param(eT(0.5)*df, eT(2)));
    }
  #else
    {
    arma_ignore(df);
    arma_ignore(n_rows);
    arma_ignore(n_cols);
    
    arma_stop_logic_error("chi2rnd(): C++11 compiler required");
    
    return obj_type();
    }
  #endif
  }



template<typename obj_type>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_Mat<obj_type>::value && is_real<typename obj_type::elem_type>::value),
  obj_type
  >::result
chi2rnd(const typename obj_type::elem_type df, const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_CXX11)
    {
    typedef typename obj_type::elem_type eT;
    
    return randg<obj_type>(s.n_rows, s.n_cols, distr_param(eT(0.5)*df, eT(2)));
    }
  #else
    {
    arma_ignore(df);
    arma_ignore(s);
    
    arma_stop_logic_error("chi2rnd(): C++11 compiler required");
    
    return obj_type();
    }
  #endif
  }



template<typename obj_type>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_Mat<obj_type>::value && is_real<typename obj_type::elem_type>::value),
  obj_type
  >::result
chi2rnd(const typename obj_type::elem_type df, const uword n_elem)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_CXX11)
    {
    typedef typename obj_type::elem_type eT;
    
    if(is_Row<obj_type>::value == true)
      {
      return randg<obj_type>(1, n_elem, distr_param(eT(0.5)*df, eT(2)));
      }
    else
      {
      return randg<obj_type>(n_elem, 1, distr_param(eT(0.5)*df, eT(2)));
      }
    }
  #else
    {
    arma_ignore(df);
    arma_ignore(n_elem);
    
    arma_stop_logic_error("chi2rnd(): C++11 compiler required");
    
    return obj_type();
    }
  #endif
  }



arma_warn_unused
inline
mat
chi2rnd(const double df, const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_CXX11)
    {
    return randg<mat>(n_rows, n_cols, distr_param(double(0.5)*df, double(2)));
    }
  #else
    {
    arma_ignore(df);
    arma_ignore(n_rows);
    arma_ignore(n_cols);
    
    arma_stop_logic_error("chi2rnd(): C++11 compiler required");
    
    return mat();
    }
  #endif
  }



arma_warn_unused
inline
mat
chi2rnd(const double df, const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_CXX11)
    {
    return randg<mat>(s.n_rows, s.n_cols, distr_param(double(0.5)*df, double(2)));
    }
  #else
    {
    arma_ignore(df);
    arma_ignore(s);
    
    arma_stop_logic_error("chi2rnd(): C++11 compiler required");
    
    return mat();
    }
  #endif
  }



arma_warn_unused
inline
vec
chi2rnd(const double df, const uword n_elem)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_CXX11)
    {
    return randg<vec>(n_elem, 1, distr_param(double(0.5)*df, double(2)));
    }
  #else
    {
    arma_ignore(df);
    arma_ignore(n_elem);
    
    arma_stop_logic_error("chi2rnd(): C++11 compiler required");
    
    return vec();
    }
  #endif
  }



//! @}
