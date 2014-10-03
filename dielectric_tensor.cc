// -------------------------------------------------------------------
// Copyright 2010 nil authors. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE NAMEPSACE EWALENA AUTHORS ``AS
// IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// NAMESPACE EWALENA AUTHORS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and
// documentation are those of the authors and should not be
// interpreted as representing official policies, either expressed or
// implied, of the namespace ewalena authors.
// -------------------------------------------------------------------

#include "dielectric_tensor.h"


namespace nil
{  

  template <enum GroupSymmetry group_symmetry, int order, typename ValueType>
  DielectricTensor<group_symmetry, order, ValueType>::DielectricTensor ()
  {}


  template <enum GroupSymmetry GroupSymm, int Order, typename ValueType>
  void 
  DielectricTensor<GroupSymm, Order, ValueType>::distribute_coefficients (std::vector<ValueType> &coefficients)
  {
    // Since template specialisations are used by this function there
    // is no need for sorting group symmetry or order of the tensor.
    // Hence, simply call the correct specialisation for the
    // non-member function in piezoelectric_tensor_base.h.
    distribute_coefficients_ (*this, coefficients);
  }
   
} // namespace nil


// -------------- Explicit Instantiations -------------------------------

// First-order tensors
template class 
nil::DielectricTensor<nil::GroupSymmetry::ZincBlende, 1, double>;

template class 
nil::DielectricTensor<nil::GroupSymmetry::Wurtzite, 1, double>;




