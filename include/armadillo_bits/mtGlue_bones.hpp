// Copyright (C) 2008-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup mtGlue
//! @{



template<typename out_eT, typename T1, typename T2, typename glue_type>
class mtGlue : public Base<out_eT, mtGlue<out_eT, T1, T2, glue_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  
  static const bool is_row = \
    (
       ((T1::is_row || T2::is_row) && is_glue_mixed_elem<glue_type>::value)
    || (T1::is_row && is_glue_mixed_times<glue_type>::value)
    || (T1::is_row && is_same_type<glue_type, glue_hist_default>::yes)
    || (T1::is_row && is_same_type<glue_type, glue_histc_default>::yes)
    );
  
  static const bool is_col = \
    ( 
       ((T1::is_col || T2::is_col) && is_glue_mixed_elem<glue_type>::value)
    || (T2::is_col && is_glue_mixed_times<glue_type>::value)
    || (T1::is_col && is_same_type<glue_type, glue_hist_default>::yes)
    || (T1::is_col && is_same_type<glue_type, glue_histc_default>::yes)
    );
  
  arma_inline  mtGlue(const T1& in_A, const T2& in_B);
  arma_inline  mtGlue(const T1& in_A, const T2& in_B, const uword in_aux_uword);
  arma_inline ~mtGlue();
  
  arma_aligned const T1&   A;         //!< first operand
  arma_aligned const T2&   B;         //!< second operand
  arma_aligned       uword aux_uword; //!< storage of auxiliary data, uword format
  };



//! @}
