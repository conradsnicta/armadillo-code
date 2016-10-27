// Copyright (C) 2008-2010 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup GlueCube
//! @{



//! analog of the Glue class, intended for Cube objects
template<typename T1, typename T2, typename glue_type>
class GlueCube : public BaseCube<typename T1::elem_type, GlueCube<T1, T2, glue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;

  arma_inline  GlueCube(const BaseCube<typename T1::elem_type, T1>& in_A, const BaseCube<typename T1::elem_type, T2>& in_B);
  arma_inline ~GlueCube();
  
  const T1& A;  //!< first operand
  const T2& B;  //!< second operand
  
  };



//! @}
