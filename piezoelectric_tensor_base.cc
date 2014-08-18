// -------------------------------------------------------------------
// @author Toby D. Young
//
// Copyright 2014 nil authors. All rights reserved.
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

#include "piezoelectric_tensor_base.h"

namespace nil
{  

  template <int rank, typename ValueType>
  PiezoelectricTensorBase<rank, ValueType>::PiezoelectricTensorBase ()
    :
    order_ (2*rank+1),
    group_symmetry_ (None)
  {}


  template <int rank, typename ValueType>
  PiezoelectricTensorBase<rank, ValueType>::PiezoelectricTensorBase (GroupSymmetry &group_symmetry)
    :
    order_ (2*rank+1),
    group_symmetry_ (group_symmetry)
  {}


  template <int rank, typename ValueType>
  PiezoelectricTensorBase<rank, ValueType>::~PiezoelectricTensorBase ()
  {}


  template <int rank, typename ValueType>
  unsigned int 
  PiezoelectricTensorBase<rank, ValueType>::order () const

  {
    return this->order_;
  }


  template <int rank, typename ValueType>
  unsigned int 
  PiezoelectricTensorBase<rank, ValueType>::dim () const
  {
    // Recall that these tensors are independent of changes in dim,
    // since they are only properly defined in 3d. Hence, return 3.
    return 3;
  }

  template <int rank, typename ValueType>
  void 
  PiezoelectricTensorBase<rank, ValueType>::reinit (GroupSymmetry &group_symmetry)
  {
    // Wipe out the tensor and re assign group symmetry.
    tensor = 0;
    this->group_symmetry_ = group_symmetry;
  }
   

} // namespace nil


// -------------- Explicit Instantiations -------------------------------

// First-order tensors
template class 
nil::PiezoelectricTensorBase<3, double>;

// second-order tensors
template class 
nil::PiezoelectricTensorBase<5, double>;

