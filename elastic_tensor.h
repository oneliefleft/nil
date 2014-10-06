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

#ifndef __nil_elastic_tensor_h
#define __nil_elastic_tensor_h

#include "group_symmetry.h"
#include "tensor_base.h"

#include <fstream>
#include <iostream>
#include <list>

namespace nil
{
  
  /**
   * \brief Elastic Tensor.
   *
   * \f$N\,\f$-order elastic tensor. 
   *
   * @author Toby D. Young  2010, 2014.
   */  
  template <enum GroupSymmetry group_symmetry, int order, typename ValueType = double>
    class ElasticTensor
    :
    public TensorBase<group_symmetry, order, 4*order, ValueType>
    {
    public:
    

    /**
     * Constructor. 
     */
    ElasticTensor ();


    /**
     * Distribute coefficients onto the matrix according to group
     * symmetry and order of the tensor.
     */
    void 
    distribute_coefficients (std::vector<ValueType> &coefficients);
    

    /* /\** */
    /*  * Read only access operator to the underlying */
    /*  * <code>deal.II</code> tensor. */
    /*  *\/     */
    /* inline */
    /* dealii::Tensor<4*order, 3, ValueType>* operator* () const */
    /* { */
    /*   return *this; */
    /* } */


    private:

    
    }; /* ElasticTensor */
  
  
  /* ----------------- Non-member functions operating on tensors. ------------ */


  /* -------------------------- Second-order tensors. ------------------------ */
  
  
  /**
   * Distribute <code>coefficients</code> on to a second-order elastic
   * tensor of zinc-blende symmetry. @note Order of the coefficients
   * is important and should be passed to this function as:
   * \f$C_{??}\f$\,.
   */  
  template <typename ValueType> 
    inline
    void 
    distribute_coefficients_ (ElasticTensor<GroupSymmetry::ZincBlende, 1, ValueType> &tensor, 
			      std::vector<ValueType>                                 &coefficients)
    {
      AssertThrow (false, dealii::ExcNotImplemented ());

      Assert (tensor.rank ()==4, dealii::ExcInternalError ());

      AssertThrow (coefficients.size ()==4,
		   dealii::ExcMessage ("The number of coefficients does not match the default number required for zinc-blende structure."));

      // Then distribute the coefficients on to the tensor. It seems
      // there is no automagic way to do this, so just insert those
      // elements that are non-zero.
      // 
      // In Voight notation these are: 
    }


  /**
   * Distribute <code>coefficients</code> on to a second-order elastic
   * tensor of wurtzite symmetry. @note Order of the coefficients is
   * important and should be passed to this function as:
   * \fC_{??}\f$\,.
   */  
  template <typename ValueType> 
    inline
    void 
    distribute_coefficients_ (ElasticTensor<GroupSymmetry::Wurtzite, 1, ValueType> &tensor, 
			      std::vector<ValueType>                               &coefficients)
    {
      Assert (tensor.rank ()==4, dealii::ExcInternalError ());

      AssertThrow (coefficients.size ()==5,
		   dealii::ExcMessage ("The number of coefficients does not match the default number required for wurtzite structure."));

      // Then distribute the coefficients on to the tensor. It seems
      // there is no automagic way to do this, so just insert those
      // elements that are non-zero.
      // 
      // In Voight notation these are: C_11 = C_22, C_12, C_13 = C_23,
      // C_33, C_44 = C_55. In total there are five independent
      // coefficients.

      // C_11 = C_22 \mapsto:
      tensor[0][0][0][0] = tensor[1][1][1][1] = coefficients[0];

      // C_12 \mapsto:
      tensor[0][0][1][1] = tensor[1][1][0][0] = coefficients[1];

      // C_13 = C_23 \mapsto:
      tensor[0][0][2][2] = tensor[1][1][2][2] = coefficients[2];

      // C_33 \mapsto:
      tensor[2][2][2][2] = coefficients[3];

      // C_44 = C55 \mapsto:
      tensor[1][2][1][2] = tensor[2][1][1][2] = tensor[2][1][2][1] = tensor[1][2][2][1] 
	=
	tensor[2][0][2][0] = tensor[0][2][2][0] = tensor[0][2][0][2] = tensor[2][0][0][2] 
	= 
	coefficients[4];

      // C_66 = C55 \mapsto:
      tensor[0][0][0][0] 
	= 
	(coefficients[0] - coefficients[1]) /2.;
      
    }


  
} /* namespace nil */


#endif /* __nil_elastic_tensor_h */


