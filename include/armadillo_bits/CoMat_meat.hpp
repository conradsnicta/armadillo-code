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


//! \addtogroup CoMat
//! @{



template<typename eT>
inline
CoMat<eT>::~CoMat()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(map_ptr)
    {
    delete map_ptr;
    }
  
  if(arma_config::debug == true)
    {
    // try to expose buggy user code that accesses deleted objects
    map_ptr = NULL;
    }
  
  arma_type_check(( is_supported_elem_type<eT>::value == false ));
  }



template<typename eT>
inline
CoMat<eT>::CoMat()
  : n_rows (0)
  , n_cols (0)
  , n_elem (0)
  , map_ptr(NULL)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold();
  }



template<typename eT>
inline
CoMat<eT>::CoMat(const uword in_n_rows, const uword in_n_cols)
  : n_rows (in_n_rows)
  , n_cols (in_n_cols)
  , n_elem (in_n_rows * in_n_cols)
  , map_ptr(NULL)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold();
  }



template<typename eT>
inline
uword
CoMat<eT>::get_n_nonzero() const
  {
  arma_extra_debug_sigprint();
  
  return uword((*map_ptr).size());
  }



template<typename eT>
inline
void
CoMat<eT>::set_size(const uword in_n_rows, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  init_warm(in_n_rows, in_n_cols);
  }



template<typename eT>
inline
void
CoMat<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  (*map_ptr).clear();
  }



template<typename eT>
inline
void
CoMat<eT>::zeros(const uword in_n_rows, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  init_warm(in_n_rows, in_n_cols);
  
  (*map_ptr).clear();
  }



template<typename eT>
arma_inline
void
CoMat<eT>::set_val_unsafe(const uword index, const eT& in_val)
  {
  arma_extra_debug_sigprint();
  
  (*map_ptr).operator[](index) = in_val;
  }



template<typename eT>
arma_inline
void
CoMat<eT>::set_val(const uword index, const eT& in_val)
  {
  arma_extra_debug_sigprint();
  
  if(in_val != eT(0))
    {
    (*map_ptr).operator[](index) = in_val;
    }
  else
    {
    (*this).erase_val(index);
    }
  }



template<typename eT>
inline
void
CoMat<eT>::erase_val(const uword index)
  {
  arma_extra_debug_sigprint();
  
  map_type& map_ref = (*map_ptr);
  
  typename map_type::iterator it     = map_ref.find(index);
  typename map_type::iterator it_end = map_ref.end();
  
  if(it != it_end)  { map_ref.erase(it); }
  }



template<typename eT>
inline
arma_warn_unused
CoMat_val<eT>
CoMat<eT>::operator[](const uword index)
  {
  return CoMat_val<eT>(*this, index);
  }



template<typename eT>
inline
arma_warn_unused
eT
CoMat<eT>::operator[](const uword index) const
  {
  map_type& map_ref = (*map_ptr);
  
  typename map_type::const_iterator it     = map_ref.find(index);
  typename map_type::const_iterator it_end = map_ref.end();
  
  return (it != it_end) ? eT((*it).second) : eT(0);
  }



template<typename eT>
inline
arma_warn_unused
CoMat_val<eT>
CoMat<eT>::operator()(const uword index)
  {
  arma_debug_check( (index >= n_elem), "CoMat::operator(): index out of bounds" );
  
  return CoMat_val<eT>(*this, index);
  }



template<typename eT>
inline
arma_warn_unused
eT
CoMat<eT>::operator()(const uword index) const
  {
  arma_debug_check( (index >= n_elem), "CoMat::operator(): index out of bounds" );
  
  map_type& map_ref = (*map_ptr);
  
  typename map_type::const_iterator it     = map_ref.find(index);
  typename map_type::const_iterator it_end = map_ref.end();
  
  return (it != it_end) ? eT((*it).second) : eT(0);
  }



template<typename eT>
inline
arma_warn_unused
CoMat_val<eT>
CoMat<eT>::at(const uword in_row, const uword in_col)
  {
  const uword index = (n_rows * in_col) + in_row;
  
  return CoMat_val<eT>(*this, index);
  }



template<typename eT>
inline
arma_warn_unused
eT
CoMat<eT>::at(const uword in_row, const uword in_col) const
  {
  const uword index = (n_rows * in_col) + in_row;
  
  map_type& map_ref = (*map_ptr);
  
  typename map_type::const_iterator it     = map_ref.find(index);
  typename map_type::const_iterator it_end = map_ref.end();
  
  return (it != it_end) ? eT((*it).second) : eT(0);
  }



template<typename eT>
inline
arma_warn_unused
CoMat_val<eT>
CoMat<eT>::operator()(const uword in_row, const uword in_col)
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "CoMat::operator(): index out of bounds" );
  
  const uword index = (n_rows * in_col) + in_row;
  
  return CoMat_val<eT>(*this, index);
  }



template<typename eT>
inline
arma_warn_unused
eT
CoMat<eT>::operator()(const uword in_row, const uword in_col) const
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "CoMat::operator(): index out of bounds" );
  
  const uword index = (n_rows * in_col) + in_row;
  
  map_type& map_ref = (*map_ptr);
  
  typename map_type::const_iterator it     = map_ref.find(index);
  typename map_type::const_iterator it_end = map_ref.end();
  
  return (it != it_end) ? eT((*it).second) : eT(0);
  }



// this function is for debugging purposes only
template<typename eT>
inline
void
CoMat<eT>::sprandu(const uword in_n_rows, const uword in_n_cols, const double density)
  {
  arma_extra_debug_sigprint();
  
  zeros(in_n_rows, in_n_cols);
  
  const uword N = uword(density * double(n_elem));
  
  const Col<eT> values(N, fill::randu);
  const eT* values_mem = values.memptr();
  
  const Col<uword> indices = linspace< Col<uword> >(0, ((n_elem > 0) ? uword(n_elem-1) : uword(0)) , N);
  const uword* indices_mem = indices.memptr();
  
  for(uword i=0; i < N; ++i)
    {
    (*this).set_val( indices_mem[i], values_mem[i] );
    }
  }



// this function is for debugging purposes only
template<typename eT>
inline
void
CoMat<eT>::print(const std::string& extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = ARMA_DEFAULT_OSTREAM.width();
    
    ARMA_DEFAULT_OSTREAM << extra_text << '\n';
    
    ARMA_DEFAULT_OSTREAM.width(orig_width);
    }
  
  const uword n_nonzero = (*map_ptr).size();
  
  const double density = (n_elem > 0) ? ((double(n_nonzero) / double(n_elem))*double(100)) : double(0);
  
  ARMA_DEFAULT_OSTREAM
    << "[matrix size: " << n_rows << 'x' << n_cols << "; n_nonzero: " << n_nonzero
    << "; density: " << density << "%]\n\n";
  
  if(n_nonzero > 0)
    {
    CoMat_const_iterator<eT> it     = (*this).begin();
    CoMat_const_iterator<eT> it_end = (*this).end();
    
    uword row = 0;
    uword col = 0;
    
    for(; it != it_end; ++it)
      {
      const eT val = (*it);
      
      if(val != eT(0))
        {
        ARMA_DEFAULT_OSTREAM << '(' << row << ", " << col << ") ";
        ARMA_DEFAULT_OSTREAM << val << '\n';
        }
      
      ++row;
      
      if(row >= n_rows)  { row = 0; col++; }
      }
    }
  
  ARMA_DEFAULT_OSTREAM.flush();
  }



template<typename eT>
inline
CoMat<eT>::init_cold()
  {
  arma_extra_debug_sigprint();
  
  // ensure that n_elem can hold the result of (n_rows * n_cols)
  
  #if (defined(ARMA_USE_CXX11) || defined(ARMA_64BIT_WORD))
    const char* error_message = "CoMat(): requested size is too large";
  #else
    const char* error_message = "CoMat(): requested size is too large; suggest to compile in C++11 mode or enable ARMA_64BIT_WORD";
  #endif
  
  arma_debug_check
    (
      (
      ( (n_rows > ARMA_MAX_UHWORD) || (n_cols > ARMA_MAX_UHWORD) )
        ? ( (double(n_rows) * double(n_cols)) > double(ARMA_MAX_UWORD) )
        : false
      ),
    error_message
    );
  
  map_ptr = new (std::nothrow) map_type;
  
  arma_check_bad_alloc( (map_ptr == NULL), "CoMat(): out of memory" );
  }



template<typename eT>
inline
void
CoMat<eT>::init_warm(const uword in_n_rows, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  if( (n_rows == in_n_rows) && (n_cols == in_n_cols))  { return; }
  
  // ensure that n_elem can hold the result of (n_rows * n_cols)
  
  #if (defined(ARMA_USE_CXX11) || defined(ARMA_64BIT_WORD))
    const char* error_message = "CoMat(): requested size is too large";
  #else
    const char* error_message = "CoMat(): requested size is too large; suggest to compile in C++11 mode or enable ARMA_64BIT_WORD";
  #endif
  
  arma_debug_check
    (
      (
      ( (in_n_rows > ARMA_MAX_UHWORD) || (in_n_cols > ARMA_MAX_UHWORD) )
        ? ( (double(in_n_rows) * double(in_n_cols)) > double(ARMA_MAX_UWORD) )
        : false
      ),
    error_message
    );
  
  const uword new_n_elem = in_n_rows * in_n_cols;
  
  access::rw(n_rows) = in_n_rows;
  access::rw(n_cols) = in_n_cols;
  access::rw(n_elem) = new_n_elem;
  
  if(new_n_elem == 0)  { (*map_ptr).clear(); }
  }






// CoMat_val



template<typename eT>
inline
CoMat_val<eT>::CoMat_val(CoMat<eT>& in_parent, const uword in_index)
  : parent(in_parent)
  , index (in_index )
  {
  }



template<typename eT>
inline
CoMat_val<eT>::operator eT() const
  {
  const CoMat<eT>& const_parent = parent;
  
  return const_parent.operator[](index);
  }



template<typename eT>
inline
void
CoMat_val<eT>::operator=(const eT in_val)
  {
  parent.set_val(index, in_val);
  }



template<typename eT>
inline
void
CoMat_val<eT>::operator+=(const eT in_val)
  {
  map_type& map_ref = *(parent.map_ptr);
  
  if(in_val != eT(0))
    {
    eT& val = map_ref.operator[](index);  // creates the element if it doesn't exist
    
    val += in_val;
    
    if(val == eT(0))  { map_ref.erase(index); }
    }
  }



template<typename eT>
inline
void
CoMat_val<eT>::operator-=(const eT in_val)
  {
  map_type& map_ref = *(parent.map_ptr);
  
  if(in_val != eT(0))
    {
    eT& val = map_ref.operator[](index);  // creates the element if it doesn't exist
    
    val -= in_val;
    
    if(val == eT(0))  { map_ref.erase(index); }
    }
  }



template<typename eT>
inline
void
CoMat_val<eT>::operator*=(const eT in_val)
  {
  map_type& map_ref = *(parent.map_ptr);
  
  typename map_type::iterator it     = map_ref.find(index);
  typename map_type::iterator it_end = map_ref.end();
  
  if(it != it_end)
    {
    if(in_val != eT(0))
      {
      eT& val = (*it).second;
      
      val *= in_val;
      
      if(val == eT(0))  { map_ref.erase(it); }
      }
    else
      {
      map_ref.erase(it);
      }
    }
  }



template<typename eT>
inline
void
CoMat_val<eT>::operator/=(const eT in_val)
  {
  arma_check( (in_val == eT(0)), "CoMat_val::operator/=(): division by zero" );
  
  map_type& map_ref = *(parent.map_ptr);
  
  typename map_type::iterator it     = map_ref.find(index);
  typename map_type::iterator it_end = map_ref.end();
  
  if(it != it_end)
    {
    if(in_val != eT(0))
      {
      eT& val = (*it).second;
      
      val /= in_val;
      
      if(val == eT(0))  { map_ref.erase(it); }
      }
    else
      {
      map_ref.erase(it);
      }
    }
  }



//! @}
