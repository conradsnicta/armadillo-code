// Copyright (C) 2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup glue_histc
//! @{


class glue_histc
   {
   public:
   
   template<typename eT>
   inline static void apply_noalias(Mat<uword>& C, const Mat<eT>& A, const Mat<eT>& B, const uword dim);
   
   template<typename T1, typename T2>
   inline static void apply(Mat<uword>& C, const mtGlue<uword,T1,T2,glue_histc>& expr);
   };



class glue_histc_default
   {
   public:
   
   template<typename T1, typename T2>
   inline static void apply(Mat<uword>& C, const mtGlue<uword,T1,T2,glue_histc_default>& expr);
   };


//! @}
