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

#ifndef __nil_piezoelectric_coefficients_h
#define __nil_piezoelectric_coefficients_h

#include "include/nil/dielectric_tensor.h"
#include "include/nil/elastic_tensor.h"
#include "include/nil/piezoelectric_tensor.h"
#include "include/nil/spontaneous_polarization_tensor.h"


namespace nil
{
  
  enum UpdateFlags
  {

    first_order_elastic,
    first_order_dielectric,
    first_order_piezoelectric,
    first_order_spontaneous_polarization,
    
    second_order_elastic,
    second_order_dielectric,
    second_order_piezoelectric,
    second_order_spontaneous_polarization
    
  }; // enum UpdateFlags
  

  /**
   * A class that handles how piezoelectric parameters are handled.
   *
   * @note Supported tensor types are: Elastic, Dielectric,
   * Piezoelectric, and SpontaneousPolarization.
   *
   * @author Toby D. Young 2014.
   */
  template <enum nil::GroupSymmetry group_symmetry, typename ValueType = double>
    class PiezoelectricCoefficients
    {
    public:
    
    
    /**
     * Constructor. This initializes an empty state with the flags
     * FirstOrder=None and SecondOrder = None.
     */
    PiezoelectricCoefficients ();
    
    
    /**
     * Constructor.
     */
    PiezoelectricCoefficients (const nil::UpdateFlags flags);
    
    
    /**
     * Destructor.
     */
    ~PiezoelectricCoefficients ();
    
    
    /**
     * Distribute coefficients onto the tensors accoring to the
     * specified group symmetry.
     */
    void distribute_coefficients ();    
    
    
    private:
    
    
    /**
     * First-order elastic tensor.
     */
    nil::ElasticTensor<group_symmetry, 1, ValueType> first_order_elastic_tensor;
    
    
    /**
     * First-order dielectric tensor.
     */
    nil::DielectricTensor<group_symmetry, 1, ValueType> first_order_dielectric_tensor;
    
    
    /**
     * First-order piezoelectric tensor.
     */
    nil::PiezoelectricTensor<group_symmetry, 1, ValueType> first_order_piezoelectric_tensor;
    
    
    /**
     * First-order spontaneous polarization tensor.
     */
    nil::SpontaneousPolarizationTensor<group_symmetry, 1, ValueType> first_order_spontaneous_polarization_tensor;
    
    
    /**
     * First-order piezoelectric tensor.
     */
    nil::PiezoelectricTensor<group_symmetry, 2, ValueType> second_order_piezoelectric_tensor;
    
    
    /**
     * First-order elastic coefficients.
     */
    std::vector<ValueType> first_order_elastic_coefficients;
    
    
    /**
     * First-order dielectric coefficients.
     */
    std::vector<ValueType> first_order_dielectric_coefficients;
    

    /**
     * First-order piezoelectric coefficients.
     */
    std::vector<ValueType> first_order_piezoelectric_coefficients;
    
    
    /**
     * First-order spontaneous polarization coefficients.
     */
    std::vector<ValueType> first_order_spontaneous_polarization_coefficients;
    
    /**
     * A local copy of the update flags handed to the constructor.
     */
    UpdateFlags update_flags;
    
    };
  
  
} // namespace nil

#endif // __nil_piezoelectric_coefficients_h
