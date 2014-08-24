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

#include <deal.II/base/tensor.h>

#include "group_symmetry.h"

#include <fstream>
#include <iostream>
#include <list>

namespace nil
{

  /**
   * \brief Piezoelectric Tensor Base.
   *
   * This is the base class for an \f$N\,\f$-order piezoelectric
   * tensor.
   *
   * @note In the inline documentation, mappings of the piezoelectric
   * constants in Voight notation is taken from:
   *
   * H. Grimmer, "The piezoelectric effect of second order in stress
   * or strain: its form for crystals and quasicrystals of any
   * symmetry", Acta Cryst. (2007) A36, 441.
   *
   * @note In the inline documentation, Various mappings from Voight
   * notation to proper tensor notation is taken from:
   *
   * J. F. Nye, "Własnosści fizyczne kryształów w Ujęciu tensorowym i
   * macierzowym" Państwowe Wydawnictwo Naukowa (1962). Tłumaczył z
   * języka angielskiego J. Rauułuszhiewicz.
   *
   * @author Toby D. Young 2014.
   */  
  template <nil::GroupSymmetry group_symmetry, int rank, typename ValueType = double>
    class PiezoelectricTensorBase
    :
    public dealii::Tensor<rank, 3, ValueType>
    {
    public:
    

    /**
     * Constructor. 
     */
    PiezoelectricTensorBase ();
    

    /**
     * Destructor. 
     */
    ~PiezoelectricTensorBase (); 


    /**
     * Distribute coefficients on to the tensor.
     */
    void distribute_coefficients ();
    

    /**
     * Reinitialise (zero out) this tensor with this symmetry.
     */
    void reinit ();


    /**
     * Return the number of non-zero elements this tensor has.
     */
    /* unsigned int n_nonzero_elements (); */


    /**
     * Make the order of this tensor public.
     */
    unsigned int order () const;
    

    /**
     * Make the dimension of this tensor public. @note This function
     * just returns the integer three, since these tensors are defined
     * in three-dimensions only.
     */
    unsigned int dim () const;


    /* /\** */
    /*  * Make the group symmetry of this tensor public. */
    /*  *\/ */
    /* std::string group_symmetry () const; */
   

    protected:
    

    /**
     * Make the order of this tensor known to all derived classes.
     */
    const int order_; 
    

    /**
     * Make the group symmetry of this tensor known to all derived classes.
     */
    GroupSymmetry group_symmetry_; 
    

    /**
     * Zero out a tensor. @note Only zero is allowed as an input to this function
     */
    /* operator = ValueType; */
    

    private:


    /**
     * The underlying tensor.
     */
    dealii::Tensor<rank, 3, ValueType> tensor; 

    }; /* PiezoelectricTensorBase */

  
  /* ----------------- Non-member functions operating on tensors. ------------ */
  
  
} /* namespace nil */


#endif /* __nil_piezoelectric_tensor_base_h */


