// Copyright (C) 2013-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup glue_min
//! @{



class glue_min
  {
  public:
  
  // dense matrices
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_min>& X);
  
  template<typename eT, typename T1, typename T2> inline static void apply(Mat< eT              >& out, const Proxy<T1>& PA, const Proxy<T2>& PB);
  
  template<typename  T, typename T1, typename T2> inline static void apply(Mat< std::complex<T> >& out, const Proxy<T1>& PA, const Proxy<T2>& PB);
  
  
  // cubes
  
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const GlueCube<T1,T2,glue_min>& X);
  
  template<typename eT, typename T1, typename T2> inline static void apply(Cube< eT              >& out, const ProxyCube<T1>& PA, const ProxyCube<T2>& PB);
  
  template<typename  T, typename T1, typename T2> inline static void apply(Cube< std::complex<T> >& out, const ProxyCube<T1>& PA, const ProxyCube<T2>& PB);
  };



//! @}

