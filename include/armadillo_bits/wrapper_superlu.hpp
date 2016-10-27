// Copyright (C) 2015 National ICT Australia (NICTA)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by Ryan Curtin



#if defined(ARMA_USE_SUPERLU)

//! \namespace superlu namespace for SuperLU functions
namespace superlu
  {
  
  template<typename eT>
  inline
  void
  gssv(superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r, SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t* stat, int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    if(is_float<eT>::value)
      {
      arma_wrapper(sgssv)(options, A, perm_c, perm_r, L, U, B, stat, info);
      }
    else
    if(is_double<eT>::value)
      {
      arma_wrapper(dgssv)(options, A, perm_c, perm_r, L, U, B, stat, info);
      }
    else
    if(is_supported_complex_float<eT>::value)
      {
      arma_wrapper(cgssv)(options, A, perm_c, perm_r, L, U, B, stat, info);
      }
    else
    if(is_supported_complex_double<eT>::value)
      {
      arma_wrapper(zgssv)(options, A, perm_c, perm_r, L, U, B, stat, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  gssvx(
        superlu_options_t* opts,
        SuperMatrix* A,
        int* perm_c, int* perm_r,
        int* etree, char* equed,
        typename get_pod_type<eT>::result* R, typename get_pod_type<eT>::result* C,
        SuperMatrix* L, SuperMatrix* U,
        void* work, int lwork,
        SuperMatrix* B, SuperMatrix* X,
        typename get_pod_type<eT>::result* rpg, typename get_pod_type<eT>::result* rcond,
        typename get_pod_type<eT>::result* ferr, typename get_pod_type<eT>::result* berr,
        GlobalLU_t* glu, mem_usage_t* mu, SuperLUStat_t* stat, int* info
       )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    if(is_float<eT>::value)
      {
      typedef float T;
      arma_wrapper(sgssvx)(opts, A, perm_c, perm_r, etree, equed, (T*)R, (T*)C, L, U, work, lwork, B, X, (T*)rpg, (T*)rcond, (T*)ferr, (T*)berr, glu, mu, stat, info);
      }
    else
    if(is_double<eT>::value)
      {
      typedef double T;
      arma_wrapper(dgssvx)(opts, A, perm_c, perm_r, etree, equed, (T*)R, (T*)C, L, U, work, lwork, B, X, (T*)rpg, (T*)rcond, (T*)ferr, (T*)berr, glu, mu, stat, info);
      }
    else
    if(is_supported_complex_float<eT>::value)
      {
      typedef float T;
      arma_wrapper(cgssvx)(opts, A, perm_c, perm_r, etree, equed, (T*)R, (T*)C, L, U, work, lwork, B, X, (T*)rpg, (T*)rcond, (T*)ferr, (T*)berr, glu, mu, stat, info);
      }
    else
    if(is_supported_complex_double<eT>::value)
      {
      typedef double T;
      arma_wrapper(zgssvx)(opts, A, perm_c, perm_r, etree, equed, (T*)R, (T*)C, L, U, work, lwork, B, X, (T*)rpg, (T*)rcond, (T*)ferr, (T*)berr, glu, mu, stat, info);
      }
    }
  
  
  
  inline
  void
  init_stat(SuperLUStat_t* stat)
    {
    arma_wrapper(StatInit)(stat);
    }


  inline
  void
  free_stat(SuperLUStat_t* stat)
    {
    arma_wrapper(StatFree)(stat);
    }
  
  
  
  inline
  void
  set_default_opts(superlu_options_t* opts)
    {
    arma_wrapper(set_default_options)(opts);
    }
  
  
  
  inline
  void
  destroy_supernode_mat(SuperMatrix* a)
    {
    arma_wrapper(Destroy_SuperNode_Matrix)(a);
    }



  inline
  void
  destroy_compcol_mat(SuperMatrix* a)
    {
    arma_wrapper(Destroy_CompCol_Matrix)(a);
    }



  inline
  void
  destroy_dense_mat(SuperMatrix* a)
    {
    arma_wrapper(Destroy_SuperMatrix_Store)(a);
    }
  
  
  
  inline
  void*
  malloc(size_t N)
    {
    return arma_wrapper(superlu_malloc)(N);
    }
  
  
  
  inline
  void
  free(void* mem)
    {
    arma_wrapper(superlu_free)(mem);
    }
  
  } // namespace superlu

#endif
