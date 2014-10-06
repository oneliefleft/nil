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

#ifndef __nil_spontaneous_polarization_tensor_h
#define __nil_spontaneous_polarization_tensor_h

#include "group_symmetry.h"
#include "tensor_base.h"

#include <fstream>
#include <iostream>
#include <list>

namespace nil
{
  
  /**
   * \brief Spontaneous polarization Tensor.
   *
   * \f$N\,\f$-order spontaneous polarization tensor. 
   *
   * @author Toby D. Young  2011, 2014.
   */  
  template <enum GroupSymmetry group_symmetry, int order, typename ValueType = double>
    class SpontaneousPolarizationTensor
    :
    public TensorBase<group_symmetry, order, order, ValueType>
    {
    public:
    

    /**
     * Constructor. 
     */
    SpontaneousPolarizationTensor ();


    /**
     * Distribute coefficients onto the matrix according to group
     * symmetry and order of the tensor.
     */
    void 
    distribute_coefficients (std::vector<ValueType> &coefficients);
    

    private:
    
    
    }; /* SpontaneousPolarizationTensor */
  
  
  /* ----------------- Non-member functions operating on tensors. ------------ */

  /* -------------------------- First-order tensors. ------------------------- */


  /**
   * Distribute <code>coefficients</code> on to a first-order
   * spontaneous polarization tensor of zinc-blende symmetry. @note
   * Zinc-blende structure has no spontaneous polarization, hence,
   * this function simply preforms sanity checks and returns.
   */
  template <typename ValueType> 
    inline
    void 
    distribute_coefficients_ (SpontaneousPolarizationTensor<GroupSymmetry::ZincBlende, 1, ValueType> &tensor, 
			      std::vector<ValueType>                                                 &coefficients)
    {
      
      Assert (tensor.rank ()==1, dealii::ExcInternalError ());
      
      AssertThrow (coefficients.size ()==0,
		   dealii::ExcMessage ("The number of coefficients does not match the default number required for zinc-blende structure."));

    }

  
  /**
   * Distribute <code>coefficients</code> on to a first-order
   * spontaneous polarization tensor of wurtzite symmetry. @note Order
   * of the coefficients is important and should be passed to this
   * function as: \f$P_3\f$\,.
   */
  template <typename ValueType> 
    inline
    void 
    distribute_coefficients_ (SpontaneousPolarizationTensor<GroupSymmetry::Wurtzite, 1, ValueType> &tensor, 
			      std::vector<ValueType>                                               &coefficients)
    {
      Assert (tensor.rank ()==1, dealii::ExcInternalError ()); 
      
      AssertThrow (coefficients.size ()==1,
		   dealii::ExcMessage ("The number of coefficients does not match the default number required for wurtzite structure."));

      // Then distribute the coefficients on to the tensor. It seems
      // there is no automagic way to do this, so just insert those
      // elements that are non-zero.
      // 
      // In Voight notation these are: P_3.

      // P_3 
      tensor[2] = coefficients[0];

    }
 
  
} /* namespace nil */


#endif /* __nil_spontaneous_polarization_tensor_h */


