// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup subview_cube_each
//! @{



template<typename eT>
class subview_cube_each_common
  {
  public:
  
  const Cube<eT>& P;
  
  inline void check_size(const Mat<eT>& A) const;
  
  
  protected:
  
  arma_inline subview_cube_each_common(const Cube<eT>& in_p);
  
  arma_cold inline const std::string incompat_size_string(const Mat<eT>& A) const;
  
  
  private:
  
  subview_cube_each_common();
  };




template<typename eT>
class subview_cube_each1 : public subview_cube_each_common<eT>
  {
  protected:
  
  arma_inline subview_cube_each1(const Cube<eT>& in_p);
  
  
  public:
  
  inline ~subview_cube_each1();
  
  // deliberately returning void
  template<typename T1> inline void operator=  (const Base<eT,T1>& x);
  template<typename T1> inline void operator+= (const Base<eT,T1>& x);
  template<typename T1> inline void operator-= (const Base<eT,T1>& x);
  template<typename T1> inline void operator%= (const Base<eT,T1>& x);
  template<typename T1> inline void operator/= (const Base<eT,T1>& x);
  template<typename T1> inline void operator*= (const Base<eT,T1>& x);
  
  
  private:
  
  friend class Cube<eT>;
  };



template<typename eT, typename TB>
class subview_cube_each2 : public subview_cube_each_common<eT>
  {
  protected:
  
  inline subview_cube_each2(const Cube<eT>& in_p, const Base<uword, TB>& in_indices);
  
  
  public:
  
  const Base<uword, TB>& base_indices;
  
  inline void check_indices(const Mat<uword>& indices) const;
  inline ~subview_cube_each2();
  
  // deliberately returning void
  template<typename T1> inline void operator=  (const Base<eT,T1>& x);
  template<typename T1> inline void operator+= (const Base<eT,T1>& x);
  template<typename T1> inline void operator-= (const Base<eT,T1>& x);
  template<typename T1> inline void operator%= (const Base<eT,T1>& x);
  template<typename T1> inline void operator/= (const Base<eT,T1>& x);
  
  
  private:
  
  friend class Cube<eT>;
  };



class subview_cube_each1_aux
  {
  public:
  
  template<typename eT, typename T2>
  static inline Cube<eT> operator_plus(const subview_cube_each1<eT>& X, const Base<eT,T2>& Y);
    
  template<typename eT, typename T2>
  static inline Cube<eT> operator_minus(const subview_cube_each1<eT>& X, const Base<eT,T2>& Y);
  
  template<typename T1, typename eT>
  static inline Cube<eT> operator_minus(const Base<eT,T1>& X, const subview_cube_each1<eT>& Y);
  
  template<typename eT, typename T2>
  static inline Cube<eT> operator_schur(const subview_cube_each1<eT>& X, const Base<eT,T2>& Y);
  
  template<typename eT, typename T2>
  static inline Cube<eT> operator_div(const subview_cube_each1<eT>& X,const Base<eT,T2>& Y);
  
  template<typename T1, typename eT>
  static inline Cube<eT> operator_div(const Base<eT,T1>& X, const subview_cube_each1<eT>& Y);
  
  template<typename eT, typename T2>
  static inline Cube<eT> operator_times(const subview_cube_each1<eT>& X,const Base<eT,T2>& Y);
  
  template<typename T1, typename eT>
  static inline Cube<eT> operator_times(const Base<eT,T1>& X, const subview_cube_each1<eT>& Y);
  };



class subview_cube_each2_aux
  {
  public:
  
  template<typename eT, typename TB, typename T2>
  static inline Cube<eT> operator_plus(const subview_cube_each2<eT,TB>& X, const Base<eT,T2>& Y);
    
  template<typename eT, typename TB, typename T2>
  static inline Cube<eT> operator_minus(const subview_cube_each2<eT,TB>& X, const Base<eT,T2>& Y);
  
  template<typename T1, typename eT, typename TB>
  static inline Cube<eT> operator_minus(const Base<eT,T1>& X, const subview_cube_each2<eT,TB>& Y);
  
  template<typename eT, typename TB, typename T2>
  static inline Cube<eT> operator_schur(const subview_cube_each2<eT,TB>& X, const Base<eT,T2>& Y);
  
  template<typename eT, typename TB, typename T2>
  static inline Cube<eT> operator_div(const subview_cube_each2<eT,TB>& X, const Base<eT,T2>& Y);
  
  template<typename T1, typename eT, typename TB>
  static inline Cube<eT> operator_div(const Base<eT,T1>& X, const subview_cube_each2<eT,TB>& Y);
  
  // TODO: operator_times
  };



//! @}
