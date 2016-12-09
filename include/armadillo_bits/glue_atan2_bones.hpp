// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup glue_atan2
//! @{



class glue_atan2
  {
  public:
  
  
  // matrices
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_atan2>& expr);
  
  template<typename T1, typename T2> inline static void apply_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P1, const Proxy<T2>& P2);
  
  
  // cubes
  
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const GlueCube<T1, T2, glue_atan2>& expr);
  
  template<typename T1, typename T2> inline static void apply_noalias(Cube<typename T1::elem_type>& out, const ProxyCube<T1>& P1, const ProxyCube<T2>& P2);
  };



//! @}
