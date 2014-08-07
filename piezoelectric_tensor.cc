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

  template <int order, typename ValueType>
  PiezoelectricTensor<order, ValueType>::PiezoelectricTensor ()
    :
    order_ (order)
  {}


  template <int order, typename ValueType>
  PiezoelectricTensor<order, ValueType>::~PiezoelectricTensor ()
  {}


  // template <int order, typename ValueType>
  // void
  // PiezoelectricTensor<order, ValueType>::distribute_coefficients (const std::vector<ValueType> &coefficients)
  // {

  //   switch (order)
  //     {
  //     case 1:
  // 	distribute_first_order_coefficients ((*this), coefficients);
  // 	break;
	
  //     case 2:
  // 	distribute_second_order_coefficients ((*this), coefficients);
  // 	break;
	
  //     default:
  // 	AssertThrow (false, dealii::ExcNotImplemented ());
  // 	break;
  //     }

  // } // PiezoelectricTensor
   
} // namespace nil


// -------------- Explicit Instantiations -------------------------------

// First-order tensors
template class 
nil::PiezoelectricTensor<1, double>;

// second-order tensors
template class 
nil::PiezoelectricTensor<2, double>;

// extern template return-type name < argument-list > ( parameter-list ) ;  (since C++11)

// extern template 
// void
// nil::distribute_coefficients<double> (nil::PiezoelectricTensor<1, double> &, 
// 				      std::vector<double>                 &);
