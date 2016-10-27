// Copyright (C) 2010-2012 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by Stanislav Funiak



//! \addtogroup span
//! @{


struct span_alt {};


template<typename Dummy = int>
class span_base
  {
  public:
  static const span_alt all;
  };


template<typename Dummy>
const span_alt span_base<Dummy>::all = span_alt();


class span : public span_base<>
  {
  public:

  uword a;
  uword b;
  bool  whole;
  
  inline
  span()
    : whole(true)
    {
    }
  
  
  inline
  span(const span_alt&)
    : whole(true)
    {
    }
  
  
  inline
  explicit
  span(const uword in_a)
    : a(in_a)
    , b(in_a)
    , whole(false)
    {
    }
  
  
  // the "explicit" keyword is required here to prevent a C++11 compiler
  // automatically converting {a,b} into an instance of span() when submatrices are specified
  inline
  explicit
  span(const uword in_a, const uword in_b)
    : a(in_a)
    , b(in_b)
    , whole(false)
    {
    }

  };



//! @}
