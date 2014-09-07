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
// implied, of the namespace ewalena authors.
// -------------------------------------------------------------------

#ifndef __nil_strain_tensor_base_h
#define __nil_strain_tensor_base_h

#include "group_symmetry.h"
#include "tensor_base.h"

#include <fstream>
#include <iostream>
#include <list>

namespace nil
{

  /**
   * \brief This is the base class for an \f$N\,\f$-order strain
   * tensor.
   *
   * @author Toby D. Young 2014.
   */  
  template <int Order, typename ValueType = double>
    class StrainTensorBase
    :
    public dealii::Tensor<2*Order, 3, ValueType>
    {
    public:
    

    /**
     * Constructor. 
     */
    StrainTensorBase ();


    /**
     * Reinitialise (zero out) this tensor.
     */
    void reinit ();
    

    private:
    
    
    /**
     * The underlying tensor.
     */
    dealii::Tensor<2*Order, 3, ValueType> tensor;
    
    }; /* StrainTensorBase */

  
  /* ----------------- Non-member functions operating on tensors. ------------ */
  
  
} /* namespace nil */


#endif /* __nil_strain_tensor_base_h */


