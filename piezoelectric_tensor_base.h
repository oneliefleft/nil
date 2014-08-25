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

#ifndef __nil_piezoelectric_tensor_base_h
#define __nil_piezoelectric_tensor_base_h

#include "group_symmetry.h"
#include "tensor_base.h"

#include <fstream>
#include <iostream>
#include <list>

namespace nil
{

  /**
   * \brief This is the base class for an \f$N\,\f$-order
   * piezoelectric tensor.
   *
   * @author Toby D. Young 2014.
   */  
  template <nil::GroupSymmetry GroupSymm, int Order, typename ValueType = double>
    class PiezoelectricTensorBase
    :
    public nil::TensorBase<GroupSymm, Order, 2*Order+1, ValueType>
    {
    public:
    

    /**
     * Constructor. 
     */
    PiezoelectricTensorBase ();
    

    private:


    /**
     * The underlying tensor.
     */
    nil::TensorBase<GroupSymm, Order, 2*Order+1, ValueType> tensor;

    }; /* PiezoelectricTensorBase */

  


  /* -------------------------- First-order tensors. ------------------------- */


  /**
   * Distribute <code>coefficients</code> on to a first-order
   * piezoelectric tensor of zinc-blende symmetry.
   */
  template <typename ValueType> 
    inline
    void 
    distribute_coefficients_ (PiezoelectricTensorBase<GroupSymmetry::ZincBlende, 1, ValueType> &tensor, 
			      std::vector<ValueType>                                           &coefficients)
    {
      
      Assert (tensor.rank ()==3, dealii::ExcInternalError ());
      
      AssertThrow (coefficients.size ()==1,
		   dealii::ExcMessage ("The number of coefficients does not match the default number required for zinc-blende structure."));
      
      // Then distribute the coefficients on to the tensor. It seems
      // there is no automagic way to do this, so just insert those
      // elements that are non-zero.
      // 
      // In Voight notation these are:  e_14 = e_26 = e_36.
      
      // e_14 \mapsto e_123 = e_132
      tensor[0][1][2] = tensor[0][2][1] = coefficients[0];
      
      // e_26 \mapsto e_212 = e_221
      tensor[1][0][1] = tensor[1][1][0] = coefficients[0];
      
      // e_36 \mapsto e_312 = e_321
      tensor[2][0][1] = tensor[2][1][0] = coefficients[0];
    }

  
  /**
   * Distribute <code>coefficients</code> on to a first-order
   * piezoelectric tensor of wurtzite symmetry.
   */
  template <typename ValueType> 
    inline
    void 
    distribute_coefficients_ (PiezoelectricTensorBase<GroupSymmetry::Wurtzite, 1, ValueType> &tensor, 
			      std::vector<ValueType>                                         &coefficients)
    {
      Assert (tensor.rank ()==3, dealii::ExcInternalError ()); 
      
      AssertThrow (coefficients.size ()==3,
		   dealii::ExcMessage ("The number of coefficients does not match the default number required for zinc-blende structure."));
    }
  
  
  /* -------------------------- Second-order tensors. ------------------------ */
  
  
  /**
   * Distribute <code>coefficients</code> on to a second-order
   * piezoelectric tensor of zinc-blende symmetry.
   */  
  template <typename ValueType> 
    inline
    void 
    distribute_coefficients_ (PiezoelectricTensorBase<GroupSymmetry::ZincBlende, 2, ValueType> &tensor, 
			      std::vector<ValueType>                                           &coefficients)
    {
      Assert (tensor.rank ()==5, dealii::ExcInternalError ()); 
      
      AssertThrow (coefficients.size ()==3,
		   dealii::ExcMessage ("The number of coefficients does not match the default number required for zinc-blende structure."));
      
      // Then distribute the coefficients on to the tensor. It seems
      // there is no automagic way to do this, so just insert those
      // elements that are non-zero.
      //
      // In Voight notation these are: e_114 = e_124 = e_156, and
      // additionally, cyclic permutations x->y->z. In total there are
      // 24 non-zero elements.
      
      // e_114 = e_225 = e_336 \mapsto:
      tensor[0][0][0][1][2] = tensor[0][0][0][2][1] = tensor[0][1][2][0][0] = tensor[0][2][1][0][0] 
	= 
	tensor[1][1][1][0][2] = tensor[1][1][1][2][0] = tensor[1][0][2][1][1] = tensor[1][2][0][1][1] 
	= 
	tensor[2][2][2][0][1] = tensor[2][2][2][1][0] = tensor[2][0][1][2][2] = tensor[2][1][0][2][2] 
	=
	coefficients[0];
      
      // e_124 = e_134 = e_215 = e_235 = e_316 = e_326 \mapsto:
      tensor[0][1][1][1][2] = tensor[0][1][1][2][1] = tensor[0][1][2][1][1] = tensor[0][2][1][1][1] 
	= 
	tensor[0][1][2][2][2] = tensor[0][2][1][2][2] = tensor[0][2][2][2][1] = tensor[0][2][2][1][2] 
	= 
	tensor[1][2][0][0][0] = tensor[1][0][2][0][0] = tensor[1][0][0][2][0] = tensor[1][0][0][0][2] 
	=
	tensor[1][2][2][2][0] = tensor[1][2][2][0][2] = tensor[1][2][0][2][2] = tensor[1][0][2][2][2] 
	=
	tensor[2][0][1][0][0] = tensor[2][1][0][0][0] = tensor[2][0][0][0][1] = tensor[2][0][0][1][0] 
	=
	tensor[2][0][1][1][1] = tensor[2][1][0][1][1] = tensor[2][1][1][0][1] = tensor[2][1][1][1][0] 
	=
	coefficients[1];
      
      // e_345 = e_246 = e_156 \mapsto:
      tensor[0][2][0][0][1] = tensor[0][0][2][0][1] = tensor[0][2][0][1][0] = tensor[0][2][0][0][1] = 
	tensor[0][0][1][2][0] = tensor[0][1][0][2][0] = tensor[0][0][1][2][0] = tensor[0][0][1][0][2] 
	= 
	tensor[1][1][2][0][1] = tensor[1][2][1][0][1] = tensor[1][1][2][1][0] = tensor[1][2][1][1][0] = 
	tensor[1][0][1][1][2] = tensor[1][1][0][1][2] = tensor[1][0][1][2][1] = tensor[1][1][0][2][1] 
	=
	tensor[2][1][2][2][0] = tensor[2][2][1][2][0] = tensor[2][1][2][0][2] = tensor[2][2][1][0][2] = 
	tensor[2][2][0][1][2] = tensor[2][0][2][1][2] = tensor[2][2][0][2][1] = tensor[2][0][2][2][1] 
	=  
	coefficients[2];
    }


  /**
   * Distribute <code>coefficients</code> on to a second-order
   * piezoelectric tensor of wurtzite symmetry.
   */  
  template <typename ValueType> 
    inline
    void 
    distribute_coefficients_ (PiezoelectricTensorBase<GroupSymmetry::Wurtzite, 2, ValueType> &tensor, 
			      std::vector<ValueType>                                         &coefficients)
    {
      Assert (tensor.rank ()==5, dealii::ExcInternalError ()); 

      AssertThrow (coefficients.size ()==8,
		   dealii::ExcMessage ("The number of coefficients does not match the default number required for zinc-blende structure."));
    }
  
} /* namespace nil */





#endif /* __nil_piezoelectric_tensor_base_h */


