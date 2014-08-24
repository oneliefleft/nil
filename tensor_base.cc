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

#include "tensor_base.h"

namespace nil
{  

  template <enum GroupSymmetry group_symmetry, int order, int rank, typename ValueType>
  TensorBase<group_symmetry, order, rank, ValueType>::TensorBase ()
    :
    group_symmetry_ (group_symmetry),
    order_          (order), 
    rank_           (rank)
  {}


  template <enum GroupSymmetry group_symmetry, int order, int rank, typename ValueType>
  TensorBase<group_symmetry, order, rank, ValueType>::~TensorBase ()
  {}


  template <enum GroupSymmetry group_symmetry, int order, int rank, typename ValueType>
  void
  TensorBase<group_symmetry, order, rank, ValueType>::distribute_coefficients ()
  {
    // virtual function called!
    AssertThrow (false, dealii::ExcInternalError ());
  }


  template <enum GroupSymmetry group_symmetry, int order, int rank, typename ValueType>
  unsigned int 
  TensorBase<group_symmetry, order, rank, ValueType>::this_order () const

  {
    return this->order_;
  }


  template <enum GroupSymmetry group_symmetry, int order, int rank, typename ValueType>
  unsigned int 
  TensorBase<group_symmetry, order, rank, ValueType>::dim () const
  {
    // Recall that these tensors are independent of changes in dim,
    // since they are only properly defined in 3d. Hence, return 3.
    return 3;
  }


  template <enum GroupSymmetry group_symmetry, int order, int rank, typename ValueType>
  void 
  TensorBase<group_symmetry, order, rank, ValueType>::reinit ()
  {
    // Wipe out the tensor 
    tensor = 0;
  } 


  template <enum GroupSymmetry group_symmetry, int order, int rank, typename ValueType>
  std::string
  TensorBase<group_symmetry, order, rank, ValueType>::this_group_symmetry () const
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
   	AssertThrow (false, dealii::ExcNotImplemented ());
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

// First-order tensors
template class 
nil::TensorBase<nil::GroupSymmetry::ZincBlende, 1, 2, double>;

// Second-order tensors
template class 
nil::TensorBase<nil::GroupSymmetry::ZincBlende, 2, 4, double>;



