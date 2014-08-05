// -------------------------------------------------------------------
// @author Toby D. Young
//
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

#ifndef __nil_piezoelectric_tensor_h
#define __nil_piezoelectric_tensor_h

#include <deal.II/base/tensor.h>

#include <fstream>
#include <iostream>
#include <list>

namespace nil
{

  /**
   * \brief Piezoelectric Tensor.
   *
   * \f$N\,\f$-order piezoelectric tensor.
   *
   * @author Toby D. Young  2010, 2011, 2014.
   */  
  template <int dim, int order, typename ValueType = double>
    class PiezoelectricTensor
    :
    public dealii::Tensor<2*order, dim, ValueType>
    {
    public:
    
    /**
     * Constructor. 
     */
    PiezoelectricTensor () 
    {}
    
    /**
     * Virtual destructor. 
     */
    virtual ~PiezoelectricTensor (); 
    
    /**
     * Distribute moduli onto this tensor field according to rules
     * governing the crystal group symmetry.
     */
    void distribute_coefficients (const std::list<ValueType> &coefficients);

    /**
     * Distribute <code>coefficients</code> on to the first-order
     * piezoelectric tensor.
     */
    void distribute_first_order_piezoelectric_coefficients (const std::list<ValueType> coefficients);
    
    /**
     * Distribute <code>coefficients</code> on to the second-order
     * piezoelectric tensor.
     */
    void distribute_second_order_piezoelectric_coefficients (const std::list<ValueType> coefficients);
    
    private:
    
    /**
     * The underlying tensor.
     */
    dealii::Tensor<2*order, dim, ValueType> tensor;
    
    }; /* PiezoelectricTensor */
  
} /* namespace nil */

#endif /* __nil_piezoelectric_tensor_h */


