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
class CoMat
  {
  public:
  
  typedef eT                                elem_type;  //!< the type of elements stored in the matrix
  typedef typename get_pod_type<eT>::result  pod_type;  //!< if eT is std::complex<T>, pod_type is T; otherwise pod_type is eT
  
  static const bool is_row = false;
  static const bool is_col = false;
  
  const uword n_rows;    //!< number of rows     (read-only)
  const uword n_cols;    //!< number of columns  (read-only)
  const uword n_elem;    //!< number of elements (read-only)
  
  
  private:
  
  typedef typename std::map<uword, eT> map_type;
  
  arma_aligned map_type* map_ptr;
  
  
  public:
  
  inline ~CoMat();
  inline  CoMat();
  inline  CoMat(const uword in_n_rows, const uword in_n_cols);
  inline  CoMat(const SizeMat& s);
  
  inline          CoMat(const CoMat<eT>& x);
  inline void operator=(const CoMat<eT>& x);
  
  inline          CoMat(const SpMat<eT>& x);
  inline void operator=(const SpMat<eT>& x);
  
  #if defined(ARMA_USE_CXX11)
  inline          CoMat(CoMat<eT>&& x);
  inline void operator=(CoMat<eT>&& x);
  #endif
  
  inline void reset();
  inline void set_size(const uword in_n_rows);
  inline void set_size(const uword in_n_rows, const uword in_n_cols);
  inline void set_size(const SizeMat& s);
  
  inline void zeros();
  inline void zeros(const uword in_n_rows);
  inline void zeros(const uword in_n_rows, const uword in_n_cols);
  inline void zeros(const SizeMat& s);
  
  inline void eye();
  inline void eye(const uword in_n_rows, const uword in_n_cols);
  inline void eye(const SizeMat& s);
  
  inline void speye();
  inline void speye(const uword in_n_rows, const uword in_n_cols);
  inline void speye(const SizeMat& s);
  
  inline arma_warn_unused CoMat_val<eT> operator[](const uword index);
  inline arma_warn_unused           eT  operator[](const uword index) const;
  
  inline arma_warn_unused CoMat_val<eT> operator()(const uword index);
  inline arma_warn_unused           eT  operator()(const uword index) const;
  
  inline arma_warn_unused CoMat_val<eT>         at(const uword in_row, const uword in_col);
  inline arma_warn_unused           eT          at(const uword in_row, const uword in_col) const;
  
  inline arma_warn_unused CoMat_val<eT> operator()(const uword in_row, const uword in_col);
  inline arma_warn_unused           eT  operator()(const uword in_row, const uword in_col) const;
  
  inline arma_warn_unused bool is_empty()  const;
  inline arma_warn_unused bool is_vec()    const;
  inline arma_warn_unused bool is_rowvec() const;
  inline arma_warn_unused bool is_colvec() const;
  inline arma_warn_unused bool is_square() const;
  
  
  inline void sprandu(const uword in_n_rows, const uword in_n_cols, const double density);
  
  inline void print(const std::string& extra_text) const;
  
  inline uword get_n_nonzero() const;
  inline void  get_locval_format(umat& locs, Col<eT>& vals) const;
  
  
  private:
  
  inline void init_cold();
  inline void init_warm(const uword in_n_rows, const uword in_n_cols);
  
  arma_inline void   set_val(const uword index, const eT& in_val);
       inline void erase_val(const uword index);
  
  inline arma_warn_unused CoMat_const_iterator<eT> begin() const;
  inline arma_warn_unused CoMat_const_iterator<eT>   end() const;
  
  friend class CoMat_val<eT>;
  friend class CoMat_const_iterator<eT>;
  };



template<typename eT>
class CoMat_val
  {
  private:
  
  arma_aligned CoMat<eT>& parent;
  
  arma_aligned const uword index;
  
  inline CoMat_val(CoMat<eT>& in_parent, const uword in_index);
  
  friend class CoMat<eT>;
  
  
  public:
  
  inline operator eT() const;
  
  inline void operator= (const eT in_val);
  inline void operator+=(const eT in_val);
  inline void operator-=(const eT in_val);
  inline void operator*=(const eT in_val);
  inline void operator/=(const eT in_val);
  
  inline void operator++();
  inline void operator++(int);
  
  inline void operator--();
  inline void operator--(int);
  };



template<typename eT>
class CoMat_const_iterator
  {
  private:
  
  arma_aligned const CoMat<eT>& parent;
  
  arma_aligned uword index;
  
  arma_aligned typename CoMat<eT>::map_type::const_iterator it;
  arma_aligned typename CoMat<eT>::map_type::const_iterator it_end;
  
  friend class CoMat<eT>;
  
  
  public:
  
  inline CoMat_const_iterator(const CoMat<eT>& in_parent, const uword in_index);
  
  inline arma_warn_unused eT operator*() const;
  
  inline void operator++();
  inline void operator++(int);
  
  inline bool operator!=(const CoMat_const_iterator<eT>& X) const;
  };



//! @}
