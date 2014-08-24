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

#include "strain_tensor_base.h"


namespace nil
{  

  template <enum GroupSymmetry group_symmetry, int rank, typename ValueType>
  StrainTensorBase<group_symmetry, rank, ValueType>::StrainTensorBase ()
    :
    order_ (rank),         /* @note this is a conversion. */
    group_symmetry_ (group_symmetry)
  {}


  template <enum GroupSymmetry group_symmetry, int rank, typename ValueType>
  StrainTensorBase<group_symmetry, rank, ValueType>::~StrainTensorBase ()
  {}


  template <enum GroupSymmetry group_symmetry, int rank, typename ValueType>
  void
  StrainTensorBase<group_symmetry, rank, ValueType>::distribute_coefficients ()
  {
    switch (group_symmetry_)
      {
      case None:
	AssertThrow (false, dealii::ExcNotImplemented ());
	break;

      case ZincBlende:
	AssertThrow (false, dealii::ExcNotImplemented ());
	break;

      case Wurtzite:
	AssertThrow (false, dealii::ExcNotImplemented ());
	break;

      default:
	AssertThrow (false, dealii::ExcNotImplemented ());
	break;
      };
  }


  template <enum GroupSymmetry group_symmetry, int rank, typename ValueType>
  unsigned int 
  StrainTensorBase<group_symmetry, rank, ValueType>::order () const

  {
    return this->order_;
  }


  template <enum GroupSymmetry group_symmetry, int rank, typename ValueType>
  unsigned int 
  StrainTensorBase<group_symmetry, rank, ValueType>::dim () const
  {
    // Recall that these tensors are independent of changes in dim,
    // since they are only properly defined in 3d. Hence, return 3.
    return 3;
  }


  template <enum GroupSymmetry group_symmetry, int rank, typename ValueType>
  void 
  StrainTensorBase<group_symmetry, rank, ValueType>::reinit ()
  {
    // Wipe out the tensor 
    tensor = 0;
  } 


  // template <enum GroupSymmetry group_symmetry, int rank, typename ValueType>
  // std::string
  // StrainTensorBase<group_symmetry, rank, ValueType>::group_symmetry () const
  // {
  //   switch (group_symmetry_)
  //   {
  //     case None:
  // 	return "None";
  // 	break;

  //     case ZincBlende:
  // 	return "ZincBlende";
  // 	break;

  //     case Wurtzite:
  // 	AssertThrow (false, dealii::ExcNotImplemented ());
  // 	break;

  //     default:
  // 	AssertThrow (false, dealii::ExcNotImplemented ());
  // 	break;
  //   };

  //   // shutup the compiler about no return value.
  //   return "";
  // } 


} // namespace nil


// -------------- Explicit Instantiations -------------------------------

// First-order tensors
template class 
nil::StrainTensorBase<nil::GroupSymmetry::ZincBlende, 2, double>;


