// -------------------------------------------------------------------
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
// implied, of the nil authors.
// -------------------------------------------------------------------

#include "../../include/nil/base/tensor_base.h"

namespace nil
{  

  template <enum GroupSymmetry GroupSymm, int Order, int Rank, typename ValueType>
  TensorBase<GroupSymm, Order, Rank, ValueType>::TensorBase ()
    :
    group_symmetry_ (GroupSymm),
    order_          (Order), 
    rank_           (Rank)
  {}


  template <enum GroupSymmetry GroupSymm, int Order, int Rank, typename ValueType>
  TensorBase<GroupSymm, Order, Rank, ValueType>::~TensorBase ()
  {}


  template <enum GroupSymmetry GroupSymm, int Order, int Rank, typename ValueType>
  unsigned int 
  TensorBase<GroupSymm, Order, Rank, ValueType>::order () const

  {
    return this->order_;
  }


  template <enum GroupSymmetry GroupSymm, int Order, int Rank, typename ValueType>
  unsigned int 
  TensorBase<GroupSymm, Order, Rank, ValueType>::rank () const
  {
    return this->rank_;
  }


  template <enum GroupSymmetry GroupSymm, int Order, int Rank, typename ValueType>
  unsigned int 
  TensorBase<GroupSymm, Order, Rank, ValueType>::dim () const
  {
    // Recall that these tensors are independent of changes in dim,
    // since they are only properly defined in 3d. Hence, return 3.
    return 3;
  }


  template <enum GroupSymmetry GroupSymm, int Order, int Rank, typename ValueType>
  void 
  TensorBase<GroupSymm, Order, Rank, ValueType>::reinit ()
  {
    // Wipe out the tensor 
    tensor = 0;
  } 


  template <enum GroupSymmetry GroupSymm, int Order, int Rank, typename ValueType>
  std::string
  TensorBase<GroupSymm, Order, Rank, ValueType>::group_symmetry () const
  {
    switch (group_symmetry_)
      {
      case None:
   	return "None";
   	break;
	
      case ZincBlende:
   	return "ZincBlende";
   	break;
	
      case Wurtzite:
   	return "Wurtzite";
   	break;
	
      default:
   	AssertThrow (false, dealii::ExcNotImplemented ());
  	break;
      };
    
    // shutup the compiler about no return value.
     return "";
  } 
    
} // namespace nil


// -------------- Explicit Instantiations -------------------------------

#include "tensor_base.inst"



