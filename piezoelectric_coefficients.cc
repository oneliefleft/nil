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

#include "piezoelectric_coefficients.h"

namespace nil
{

  template <enum nil::GroupSymmetry group_symmetry, typename ValueType>
  PiezoelectricCoefficients<group_symmetry, ValueType>::PiezoelectricCoefficients ()
  {}


  template <enum nil::GroupSymmetry group_symmetry, typename ValueType>
  PiezoelectricCoefficients<group_symmetry, ValueType>::PiezoelectricCoefficients (const nil::UpdateFlags flags)
  {
    if (flags == nil::UpdateFlags::first_order_elastic)
      {
	// FirstOrder::Elastic = true;
	// FirstOrder::None    = false;
      }


    // Piezolectric coupling without an elastic and electric field
    // make little (or no) physical sense. Do not allow it.
    // if ((FirstOrder::Piezoelectric) 
    // 	&& 
    // 	((!FirstOrder::Elastic) || (!FirstOrder::Dielectric)))
    //   AssertThrow (false, dealii::ExcInternalError ());
  }


  template <enum nil::GroupSymmetry group_symmetry, typename ValueType>
  PiezoelectricCoefficients<group_symmetry, ValueType>::~PiezoelectricCoefficients ()
  {}


  // @note, this is the part where input can be read from the
  // parameter file, right?
  template <enum nil::GroupSymmetry group_symmetry, typename ValueType>
  void
  PiezoelectricCoefficients<group_symmetry, ValueType>::distribute_coefficients ()
  {
    // if (this->first_order_elastic)
    //   first_order_elastic_tensor.distribute_coefficients (first_order_elastic_coefficients);
  }



} // namespace nil

// -------------- Explicit Instantiations -------------------------------

template class 
nil::PiezoelectricCoefficients<nil::GroupSymmetry::Wurtzite, double>;
