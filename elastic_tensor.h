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
#include "elastic_tensor_base.h"

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
    public ElasticTensorBase<group_symmetry, order, ValueType>
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
    

    private:
    
    /**
     * The underlying tensor.
     */
    ElasticTensorBase<group_symmetry, order, ValueType> tensor;
    
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
    distribute_coefficients_ (ElasticTensor<GroupSymmetry::ZincBlende, 2, ValueType> &tensor, 
			      std::vector<ValueType>                                 &coefficients)
    {
      AssertThrow (false, dealii::ExcNotImplemented ());
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
    distribute_coefficients_ (ElasticTensor<GroupSymmetry::Wurtzite, 2, ValueType> &tensor, 
			      std::vector<ValueType>                               &coefficients)
    {
      AssertThrow (false, dealii::ExcNotImplemented ());
    }


  
} /* namespace nil */


#endif /* __nil_piezoelectric_tensor_h */


