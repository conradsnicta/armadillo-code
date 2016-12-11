// Copyright (C) 2008-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup Glue
//! @{



//! Class for storing data required for delayed binary operations,
//! such as the operands (e.g. two matrices) and the binary operator (e.g. addition).
//! The operands are stored as references (which can be optimised away),
//! while the operator is "stored" through the template definition (glue_type).
//! The operands can be 'Mat', 'Row', 'Col', 'Op', and 'Glue'.
//! Note that as 'Glue' can be one of the operands, more than two matrices can be stored.
//!
//! For example, we could have: Glue<Mat, Mat, glue_times>
//! 
//! Another example is: Glue< Op<Mat, op_htrans>, Op<Mat, op_inv>, glue_times >



template<typename T1, typename T2, typename glue_type>
class Glue : public Base<typename T1::elem_type, Glue<T1, T2, glue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool is_row = \
       (is_same_type<glue_type,glue_times>::value && T1::is_row)
    || (is_same_type<glue_type,glue_conv>::value  && T1::is_row)
    || (is_same_type<glue_type,glue_join_rows>::value && T1::is_row && T2::is_row)
    || (is_same_type<glue_type,glue_atan2>::value && (T1::is_row || T2::is_row))
    || (is_same_type<glue_type,glue_hypot>::value && (T1::is_row || T2::is_row));
    
  static const bool is_col = \
       (is_same_type<glue_type,glue_times>::value && T2::is_col)
    || (is_same_type<glue_type,glue_conv>::value  && T1::is_col)
    || (is_same_type<glue_type,glue_join_cols>::value && T1::is_col && T2::is_col)
    || (is_same_type<glue_type,glue_atan2>::value && (T1::is_col || T2::is_col))
    || (is_same_type<glue_type,glue_hypot>::value && (T1::is_col || T2::is_col));
  
  arma_inline  Glue(const T1& in_A, const T2& in_B);
  arma_inline  Glue(const T1& in_A, const T2& in_B, const uword in_aux_uword);
  arma_inline ~Glue();
  
  const T1&   A;          //!< first operand
  const T2&   B;          //!< second operand
        uword aux_uword;  //!< storage of auxiliary data, uword format
  };



//! @}
