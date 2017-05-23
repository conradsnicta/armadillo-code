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


//! \addtogroup mp_misc
//! @{




template<typename eT, const bool use_smaller_thresh = false>
struct mp_len
  {
  arma_inline
  static
  bool
  test(const uword n_elem)
    {
    return (is_cx<eT>::yes || use_smaller_thresh) ? (n_elem >= (arma_config::mp_threshold/uword(2))) : (n_elem >= arma_config::mp_threshold);
    }
  };



//! @}
