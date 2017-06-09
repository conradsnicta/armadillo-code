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


//! \addtogroup mp_misc
//! @{




template<typename eT>
struct mp_allow
  {
  arma_inline
  static
  bool
  eval(const uword n_elem, const bool heavy, const bool heavy_dual = false)
    {
    #if defined(ARMA_USE_OPENMP)
      {
      const bool length_ok = (heavy)
                             ? ((is_cx<eT>::yes || heavy_dual) ? (n_elem >= (arma_config::mp_thresh_b/uword(2))) : (n_elem >= arma_config::mp_thresh_b))
                             : ((is_cx<eT>::yes              ) ? (n_elem >= (arma_config::mp_thresh_a/uword(2))) : (n_elem >= arma_config::mp_thresh_a));
      
      if(length_ok)
        {
        if(omp_in_parallel())  { return false; }
        }
      
      return length_ok;
      }
    #else
      arma_ignore(n_elem);
      arma_ignore(heavy);
      arma_ignore(heavy_dual);
      return false;
    #endif
    }
  };



struct mp_thread_limit
  {
  arma_inline
  static
  int
  get(const bool heavy)
    {
    #if defined(ARMA_USE_OPENMP)
      int n_wanted  = int( (heavy) ? int(arma_config::mp_threads) : int((std::min)(int(2), int(arma_config::mp_threads))) );
      int n_threads = int( (std::min)(int(n_wanted), int((std::max)(int(1), int(omp_get_max_threads())))) );
    #else
      arma_ignore(heavy);
      int n_threads = int(1);
    #endif
    
    return n_threads;
    }
  };



//! @}
