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



//! \addtogroup fn_repelem
//! @{


template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_repelem>
  >::result
repelem(const T1& A, const uword r, const uword c)
  {
  arma_extra_debug_sigprint();

  return Op<T1, op_repelem>(A, r, c);
  }



template<typename T1>
arma_warn_unused
arma_inline
const SpOp<T1, spop_repelem>
repelem(const SpBase<typename T1::elem_type,T1>& A, const uword r, const uword c)
  {
  arma_extra_debug_sigprint();

  return SpOp<T1, spop_repelem>(A.get_ref(), r, c);
  }



//! @}
