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


//! \addtogroup SpMat
//! @{


///////////////////////////////////////////////////////////////////////////////
// SpMat::iterator_base implementation                                       //
///////////////////////////////////////////////////////////////////////////////


template<typename eT>
inline
SpMat<eT>::iterator_base::iterator_base()
  : M(NULL)
  , internal_col(0)
  , internal_pos(0)
  {
  // Technically this iterator is invalid (it does not point to a valid element)
  }



template<typename eT>
inline
SpMat<eT>::iterator_base::iterator_base(const SpMat<eT>& in_M)
  : M(&in_M)
  , internal_col(0)
  , internal_pos(0)
  {
  // Technically this iterator is invalid (it may not point to a valid element)
  }



template<typename eT>
inline
SpMat<eT>::iterator_base::iterator_base(const SpMat<eT>& in_M, const uword in_col, const uword in_pos)
  : M(&in_M)
  , internal_col(in_col)
  , internal_pos(in_pos)
  {
  // Nothing to do.
  }



template<typename eT>
arma_inline
eT
SpMat<eT>::iterator_base::operator*() const
  {
  return M->values[internal_pos];
  }



///////////////////////////////////////////////////////////////////////////////
// SpMat::const_iterator implementation                                      //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator()
  : iterator_base()
  {
  }

template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator(const SpMat<eT>& in_M, uword initial_pos)
  : iterator_base(in_M, 0, initial_pos)
  {
  // Corner case for empty matrices.
  if(in_M.n_nonzero == 0)
    {
    iterator_base::internal_col = in_M.n_cols;
    return;
    }

  // Determine which column we should be in.
  while(iterator_base::M->col_ptrs[iterator_base::internal_col + 1] <= iterator_base::internal_pos)
    {
    iterator_base::internal_col++;
    }
  }



template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator(const SpMat<eT>& in_M, uword in_row, uword in_col)
  : iterator_base(in_M, in_col, 0)
  {
  // So we have a position we want to be right after.  Skip to the column.
  iterator_base::internal_pos = iterator_base::M->col_ptrs[iterator_base::internal_col];

  // Now we have to make sure that is the right column.
  while(iterator_base::M->col_ptrs[iterator_base::internal_col + 1] <= iterator_base::internal_pos)
    {
    iterator_base::internal_col++;
    }

  // Do we have to search for an element?  We only do this if we are in the
  // right column.
  if (iterator_base::internal_col == in_col)
    {
    const uword      col_offset = iterator_base::M->col_ptrs[iterator_base::internal_col];
    const uword next_col_offset = iterator_base::M->col_ptrs[iterator_base::internal_col + 1];

    const uword* start_ptr = &iterator_base::M->row_indices[     col_offset];
    const uword* end_ptr   = &iterator_base::M->row_indices[next_col_offset];

    // Perform a binary search.
    const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, in_row);

    const uword offset = pos_ptr - start_ptr;

    // We have to increment the column if the element was the first of the next
    // column.
    if (pos_ptr == end_ptr)
      ++iterator_base::internal_col;
    iterator_base::internal_pos = col_offset + offset;
    }
  }



template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator(const SpMat<eT>& in_M, const uword /* in_row */, const uword in_col, const uword in_pos)
  : iterator_base(in_M, in_col, in_pos)
  {
  // Nothing to do.
  }



template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator(const typename SpMat<eT>::const_iterator& other)
  : iterator_base(*other.M, other.internal_col, other.internal_pos)
  {
  // Nothing to do.
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::const_iterator&
SpMat<eT>::const_iterator::operator++()
  {
  ++iterator_base::internal_pos;

  if (iterator_base::internal_pos == iterator_base::M->n_nonzero)
    {
    iterator_base::internal_col = iterator_base::M->n_cols;
    return *this;
    }

  // Check to see if we moved a column.
  while (iterator_base::M->col_ptrs[iterator_base::internal_col + 1] <= iterator_base::internal_pos)
    {
    ++iterator_base::internal_col;
    }

  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::const_iterator
SpMat<eT>::const_iterator::operator++(int)
  {
  typename SpMat<eT>::const_iterator tmp(*this);

  ++(*this);

  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::const_iterator&
SpMat<eT>::const_iterator::operator--()
  {
  --iterator_base::internal_pos;

  // First, see if we moved back a column.
  while (iterator_base::internal_pos < iterator_base::M->col_ptrs[iterator_base::internal_col])
    {
    --iterator_base::internal_col;
    }

  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::const_iterator
SpMat<eT>::const_iterator::operator--(int)
  {
  typename SpMat<eT>::const_iterator tmp(*this);

  --(*this);

  return tmp;
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator==(const const_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator!=(const const_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator==(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator!=(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator==(const const_row_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator!=(const const_row_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator==(const row_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator!=(const row_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator==(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator!=(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



///////////////////////////////////////////////////////////////////////////////
// SpMat::iterator implementation                                            //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
arma_hot
SpValProxy<SpMat<eT> >
SpMat<eT>::iterator::operator*()
  {
  return SpValProxy<SpMat<eT> >(
    iterator_base::M->row_indices[iterator_base::internal_pos],
    iterator_base::internal_col,
    access::rw(*iterator_base::M),
    &access::rw(iterator_base::M->values[iterator_base::internal_pos]));
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::iterator&
SpMat<eT>::iterator::operator++()
  {
  const_iterator::operator++();
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::iterator
SpMat<eT>::iterator::operator++(int)
  {
  typename SpMat<eT>::iterator tmp(*this);

  const_iterator::operator++();

  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::iterator&
SpMat<eT>::iterator::operator--()
  {
  const_iterator::operator--();
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::iterator
SpMat<eT>::iterator::operator--(int)
  {
  typename SpMat<eT>::iterator tmp(*this);

  const_iterator::operator--();

  return tmp;
  }



///////////////////////////////////////////////////////////////////////////////
// SpMat::const_row_iterator implementation                                  //
///////////////////////////////////////////////////////////////////////////////

/**
 * Initialize the const_row_iterator.
 */

template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator() { }



template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator(const SpMat<eT>& in_M, uword initial_pos)
  : orig_m(&in_M)
  , trans_m(in_M.st())
  , it(trans_m, initial_pos)
  { /* Nothing to do. */ }



template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator(const SpMat<eT>& in_M, uword in_row, uword in_col)
  : orig_m(&in_M)
  , trans_m(in_M.st())
  , it(trans_m, in_col, in_row)
  { }



/**
 * Initialize the const_row_iterator from another const_row_iterator.
 */
template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator(const typename SpMat<eT>::const_row_iterator& other)
  : orig_m(other.orig_m)
  , trans_m(other.trans_m)
  {
  it.M = &trans_m;
  it.internal_pos = other.it.internal_pos;
  it.internal_col = other.it.internal_col;
  }



template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator(typename SpMat<eT>::const_row_iterator&& other)
  : orig_m(other.orig_m)
  , trans_m(std::move(other.trans_m))
  , it(std::move(other.it))
  {
  other.it.M = &trans_m;
  }



template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator(const typename SpMat<eT>::row_iterator& other)
  : orig_m(other.orig_m)
  , trans_m(other.trans_m)
  {
  it.M = &trans_m;
  it.internal_pos = other.it.internal_pos;
  it.internal_col = other.it.internal_col;
  }



template<typename eT>
inline
typename SpMat<eT>::const_row_iterator&
SpMat<eT>::const_row_iterator::operator=(const typename SpMat<eT>::const_row_iterator& other)
  {
  if (this == &other)
    return *this;

  orig_m = other.orig_m;
  trans_m = other.trans_m;
  it = other.it;
  it.M = &trans_m;

  return *this;
  }



template<typename eT>
inline
typename SpMat<eT>::const_row_iterator&
SpMat<eT>::const_row_iterator::operator=(typename SpMat<eT>::const_row_iterator&& other)
  {
  if (this == &other)
    return *this;

  orig_m = other.orig_m;
  trans_m = std::move(other.trans_m);
  it = std::move(other.it);
  it.M = &trans_m;

  return *this;
  }



/**
 * Increment the row_iterator.
 */
template<typename eT>
inline
arma_hot
typename SpMat<eT>::const_row_iterator&
SpMat<eT>::const_row_iterator::operator++()
  {
  ++it;
  return *this;
  }



/**
 * Increment the row_iterator (but do not return anything).
 */
template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::const_row_iterator
SpMat<eT>::const_row_iterator::operator++(int)
  {
  typename SpMat<eT>::const_row_iterator tmp(*this);

  ++(*this);

  return tmp;
  }



/**
 * Decrement the row_iterator.
 */
template<typename eT>
inline
arma_hot
typename SpMat<eT>::const_row_iterator&
SpMat<eT>::const_row_iterator::operator--()
  {
  --it;
  return *this;
  }



/**
 * Decrement the row_iterator.
 */
template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::const_row_iterator
SpMat<eT>::const_row_iterator::operator--(int)
  {
  typename SpMat<eT>::const_row_iterator tmp(*this);

  --(*this);

  return tmp;
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator==(const const_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator!=(const const_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator==(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator!=(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator==(const const_row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator!=(const const_row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator==(const row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator!=(const row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator==(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator!=(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



///////////////////////////////////////////////////////////////////////////////
// SpMat::row_iterator implementation                                        //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
SpMat<eT>::row_iterator::row_iterator() { }



template<typename eT>
inline
SpMat<eT>::row_iterator::row_iterator(SpMat<eT>& in_M, uword initial_pos)
  : orig_m(&in_M)
  , trans_m(in_M.st())
  , it(trans_m, initial_pos)
  { /* Nothing to do. */ }



template<typename eT>
inline
SpMat<eT>::row_iterator::row_iterator(SpMat<eT>& in_M, uword in_row, uword in_col)
  : orig_m(&in_M)
  , trans_m(in_M.st())
  , it(trans_m, in_col, in_row)
  { }



/**
 * Initialize the const_row_iterator from another const_row_iterator.
 */
template<typename eT>
inline
SpMat<eT>::row_iterator::row_iterator(const typename SpMat<eT>::row_iterator& other)
  : orig_m(other.orig_m)
  , trans_m(other.trans_m)
  {
  it.M = &trans_m;
  it.internal_pos = other.it.internal_pos;
  it.internal_col = other.it.internal_col;
  }



template<typename eT>
inline
SpMat<eT>::row_iterator::row_iterator(typename SpMat<eT>::row_iterator&& other)
  : orig_m(other.orig_m)
  , trans_m(std::move(other.trans_m))
  , it(std::move(other.it))
  {
  it.M = &trans_m;
  }



template<typename eT>
inline
typename SpMat<eT>::row_iterator&
SpMat<eT>::row_iterator::operator=(const typename SpMat<eT>::row_iterator& other)
  {
  if (this == &other)
    return *this;

  orig_m = other.orig_m;
  trans_m = other.trans_m;
  it = other.it;
  it.M = &trans_m;

  return *this;
  }



template<typename eT>
inline
typename SpMat<eT>::row_iterator&
SpMat<eT>::row_iterator::operator=(typename SpMat<eT>::row_iterator&& other)
  {
  if (this == &other)
    return *this;

  orig_m = other.orig_m;
  trans_m = std::move(other.trans_m);
  it = std::move(other.it);
  it.M = &trans_m;

  return *this;
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::row_iterator&
SpMat<eT>::row_iterator::operator++()
  {
  ++it;
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::row_iterator
SpMat<eT>::row_iterator::operator++(int)
  {
  typename SpMat<eT>::row_iterator tmp(*this);

  ++it;

  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::row_iterator&
SpMat<eT>::row_iterator::operator--()
  {
  --it;
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::row_iterator
SpMat<eT>::row_iterator::operator--(int)
  {
  typename SpMat<eT>::row_iterator tmp(*this);

  --it;

  return tmp;
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator==(const const_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator!=(const const_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator==(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator!=(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator==(const const_row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator!=(const const_row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator==(const row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator!=(const row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator==(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == col());
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::row_iterator::operator!=(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != col());
  }



//! @}
