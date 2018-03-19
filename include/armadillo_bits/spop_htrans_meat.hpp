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


//! \addtogroup spop_htrans
//! @{



template<typename T1>
arma_hot
inline
void
spop_htrans::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_htrans>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  spop_strans::apply(out, in);
  }



template<typename T1>
arma_hot
inline
void
spop_htrans::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_htrans>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type  eT;
  
  const SpProxy<T1> p(in.m);
  
  const uword N = p.get_n_nonzero();
  
  if(N == uword(0))
    {
    out.zeros(p.get_n_cols(), p.get_n_rows());
    return;
    }
  
  umat    locs(2, N);
  Col<eT> vals(   N);
  
  uword* locs_mem = locs.memptr();
  eT*    vals_mem = vals.memptr();
  
  typename SpProxy<T1>::const_iterator_type it = p.begin();
  
  for(uword i=0; i < N; ++i)
    {
    const uword row = it.col();
    const uword col = it.row();
    
    (*locs_mem) = row;  locs_mem++;
    (*locs_mem) = col;  locs_mem++;
    
    (*vals_mem) = std::conj(*it);  vals_mem++;
    
    ++it;
    }
  
  SpMat<eT> tmp(locs, vals, p.get_n_cols(), p.get_n_rows(), true, false);
  
  out.steal_mem(tmp);
  }



//! @}
