// Copyright (C) 2009-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_solve
//! @{



//
// solve_gen


template<typename T1, typename T2>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve_gen> >::result
solve
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_solve_gen>(A.get_ref(), B.get_ref(), opts.flags);
  }



template<typename T1, typename T2>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve_gen> >::result
solve
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("solve(A,B,bool) is deprecated and will be removed; change to solve(A,B)");
  
  return Glue<T1, T2, glue_solve_gen>(A.get_ref(), B.get_ref(), solve_opts::flag_none);
  }



template<typename T1, typename T2>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve_gen> >::result
solve
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("solve(A,B,char*) is deprecated and will be removed; change to solve(A,B)");
  
  return Glue<T1, T2, glue_solve_gen>(A.get_ref(), B.get_ref(), solve_opts::flag_none);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  return glue_solve_gen::apply(out, A.get_ref(), B.get_ref(), opts.flags);
  }



template<typename T1, typename T2>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("solve(X,A,B,bool) is deprecated and will be removed; change to solve(X,A,B)");
  
  return glue_solve_gen::apply(out, A.get_ref(), B.get_ref(), solve_opts::flag_none);
  }



template<typename T1, typename T2>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("solve(X,A,B,char*) is deprecated and will be removed; change to solve(X,A,B)");
  
  return glue_solve_gen::apply(out, A.get_ref(), B.get_ref(), solve_opts::flag_none);
  }



//
// solve_tri


template<typename T1, typename T2>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve_tri> >::result
solve
  (
  const Op<T1, op_trimat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = opts.flags;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return Glue<T1, T2, glue_solve_tri>(A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve_tri> >::result
solve
  (
  const Op<T1, op_trimat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("solve(A,B,bool) is deprecated and will be removed; change to solve(A,B)");
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return Glue<T1, T2, glue_solve_tri>(A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve_tri> >::result
solve
  (
  const Op<T1, op_trimat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("solve(A,B,char*) is deprecated and will be removed; change to solve(A,B)");
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return Glue<T1, T2, glue_solve_tri>(A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const   Op<T1, op_trimat>&             A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = opts.flags;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return glue_solve_tri::apply(out, A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const   Op<T1, op_trimat>&             A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("solve(X,A,B,bool) is deprecated and will be removed; change to solve(X,A,B)");
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return glue_solve_tri::apply(out, A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const   Op<T1, op_trimat>&             A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("solve(X,A,B,char*) is deprecated and will be removed; change to solve(X,A,B)");
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return glue_solve_tri::apply(out, A.m, B.get_ref(), flags);
  }



//! @}
