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

#include "piezoelectric_tensor.h"


namespace nil
{  

  template <int dim, int order, typename ValueType>
  PiezoelectricTensor<dim, order, ValueType>::~PiezoelectricTensor ()
  {}


  template <int dim, int order, typename ValueType>
  void 
  PiezoelectricTensor<dim, order, ValueType>::distribute_first_order_piezoelectric_coefficients 
  (const std::vector<ValueType> coefficients)
  {
    Assert (tensor.rank==3,
	    dealii::ExcMessage ("The rank of this tesor is not correct."));

    Assert (coefficients.size ()!=0, 
	    dealii::ExcMessage ("The number of coefficients can not be zero."));
    
    // At this point we are interested in zinc-blende structure only,
    // hence the number of independent coefficients is one.
    Assert (coefficients.size ()==1,
     	    dealii::ExcMessage ("The number of coefficients does not match the default number required for zinc-blende structure."));
    
    // Then distribute the coefficients on to the tensor. It seems
    // there is no automagic way to do this, so first zero out all the
    // elemtns and second insert those elements that are non-zero.
    //
    // In Voight notation these are:  e_14=e_26=e_36.
    
    
  }
  
  
  template <int dim, int order, typename ValueType>
  void 
  PiezoelectricTensor<dim, order, ValueType>::distribute_second_order_piezoelectric_coefficients 
  (const std::vector<ValueType> coefficients)
  {
    Assert (tensor.rank==5,
	    dealii::ExcMessage ("The rank of this tesor is not correct."));

    Assert (coefficients.size ()!=0, 
	    dealii::ExcMessage ("The number of coefficients can not be zero."));
    
    // At this point we are interested in zinc-blende structure only,
    // hence the number of independent coefficients is one.
    Assert (coefficients.size ()==3,
     	    dealii::ExcMessage ("The number of coefficients does not match the default number required for zinc-blende structure."));

    // At this point we are interested in zinc-blende structure only,
    // hence:
    Assert (coefficients.size ()==3, 
     	    dealii::ExcMessage ("The number of coefficients does not match the default number required for zinc-blende structure."));
}

  
} // namespace nil


// -------------- Explicit Instantiations -------------------------------

// First-order tensors
template class 
nil::PiezoelectricTensor<3, 1, double>;

// second-order tensors
template class 
nil::PiezoelectricTensor<3, 2, double>;

