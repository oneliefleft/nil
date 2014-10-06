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

#ifndef __nil_tensor_base_h
#define __nil_tensor_base_h

#include <deal.II/base/tensor.h>

#include "../../include/nil/group_symmetry.h"

#include <fstream>
#include <iostream>
#include <list>

namespace nil
{

  /**
   * \brief This is the base class for an \f$N\,\f$-order tensors.
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
  template <enum GroupSymmetry, int Order, int Rank, typename ValueType = double>
    class TensorBase
    :
    public dealii::Tensor<Rank, 3, ValueType>
    {
    public:
    

    /**
     * Constructor. 
     */
    TensorBase ();
    

    /**
     * Virtual destructor. 
     */
    virtual
    ~TensorBase (); 
    

    /**
     * Reinitialise (zero out) this tensor.
     */
    void reinit ();


    /**
     * Return the order of this tensor.
     */
    unsigned int order () const;


    /**
     * Return the rank of this tensor.
     */
    unsigned int rank () const;
    

    /**
     * Return the dimension of this tensor. @note This function just
     * returns the integer three, since these tensors are properly
     * defined in three-dimensions only.
     */
    unsigned int dim () const;


    /**
     * Make the group symmetry of this tensor public. 
     */
    std::string group_symmetry () const; 
   

    protected:
    

    /**
     * Make the group symmetry of this tensor known to all derived classes.
     */
    const GroupSymmetry group_symmetry_; 


    /**
     * Make the order of this tensor known to all derived classes.
     */
    const int order_; 
    

    /**
     * Make the rank of this tensor known to all derived classes.
     */
    const int rank_; 


    protected:


    /**
     * The underlying tensor.
     */
    dealii::Tensor<Rank, 3, ValueType> tensor; 

    }; /* TensorBase */

  
  /* ----------------- Non-member functions operating on tensors. ------------ */
  
  
} /* namespace nil */


#endif /* __nil_tensor_base_h */


