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


//! \addtogroup op_flip
//! @{



template<typename T1>
inline
void
op_flipud::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_flipud>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> U(in.m);
  
  op_flipud::apply_direct(out, U.M);
  }



template<typename eT>
inline
void
op_flipud::apply_direct(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword X_n_rows_m1 = X_n_rows - 1;
  
  if(&out != &X)
    {
    out.copy_size(X);
    
    if(X_n_cols == 1)
      {
      const eT*   X_mem =   X.memptr();
            eT* out_mem = out.memptr();
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        out_mem[row] = X_mem[X_n_rows_m1 - row];
        }
      }
    else
      {
      for(uword col=0; col < X_n_cols; ++col)
        {
        const eT*   X_colmem =   X.colptr(col);
              eT* out_colmem = out.colptr(col);
        
        for(uword row=0; row < X_n_rows; ++row)
          {
          out_colmem[row] = X_colmem[X_n_rows_m1 - row];
          }
        }
      }
    }
  else  // in-place operation
    {
    const uword N = X_n_rows / 2;
    
    if(X_n_cols == 1)
      {
      eT* out_mem = out.memptr();
      
      for(uword row=0; row < N; ++row)
        {
        std::swap(out_mem[row], out_mem[X_n_rows_m1 - row]);
        }
      }
    else
      {
      for(uword col=0; col < X_n_cols; ++col)
        {
        eT* out_colmem = out.colptr(col);
        
        for(uword row=0; row < N; ++row)
          {
          std::swap(out_colmem[row], out_colmem[X_n_rows_m1 - row]);
          }
        }
      }
    }
  }



template<typename T1>
inline
void
op_fliplr::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_fliplr>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> U(in.m);
  
  op_fliplr::apply_direct(out, U.M);
  }



template<typename eT>
inline
void
op_fliplr::apply_direct(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword X_n_cols = X.n_cols;
  const uword X_n_rows = X.n_rows;
  
  const uword X_n_cols_m1 = X_n_cols - 1;
  
  if(&out != &X)
    {
    out.copy_size(X);
    
    if(X_n_rows == 1)
      {
      const eT*   X_mem =   X.memptr();
            eT* out_mem = out.memptr();
      
      for(uword i=0; i < X_n_cols; ++i)  { out_mem[i] = X_mem[X_n_cols_m1 - i]; }
      }
    else
      {
      for(uword i=0; i < X_n_cols; ++i)  { out.col(i) = X.col(X_n_cols_m1 - i); }
      }
    }
  else  // in-place operation
    {
    const uword N = X_n_cols / 2;
    
    if(X_n_rows == 1)
      {
      eT* out_mem = out.memptr();
      
      for(uword i=0; i < N; ++i)  { std::swap(out_mem[i], out_mem[X_n_cols_m1 - i]); }
      }
    else
      {
      for(uword i=0; i < N; ++i)  { out.swap_cols(i, X_n_cols_m1 - i); }
      }
    }
  }



//! @}
