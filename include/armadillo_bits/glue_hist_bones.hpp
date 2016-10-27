// Copyright (C) 2012-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup glue_hist
//! @{


class glue_hist
   {
   public:
   
   template<typename eT>
   inline static void apply_noalias(Mat<uword>& out, const Mat<eT>& X, const Mat<eT>& C, const uword dim);
   
   template<typename T1, typename T2>
   inline static void apply(Mat<uword>& out, const mtGlue<uword,T1,T2,glue_hist>& expr);
   };



class glue_hist_default
   {
   public:
   
   template<typename T1, typename T2>
   inline static void apply(Mat<uword>& out, const mtGlue<uword,T1,T2,glue_hist_default>& expr);
   };


//! @}
