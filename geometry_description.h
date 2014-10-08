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

#ifndef __nil_geometry_description_h
#define __nil_geometry_description_h

#include <deal.II/base/function.h>

namespace nil
{
  
  namespace GeometryDescription
  {
    
    /**
     * A class that produces a solid hyper-cube from a scalar
     * function.  The hyper-cube volume is the tensor product interval
     * \f$[left, right]^{\text{dim}}\f$\,. By default a unit cube
     * centred on the origin is created.
     *
     * @author Toby D. Young 2011, 2014
     */
    template <int dim, typename ValueType = double>
      class HyperCube
      :
      public dealii::Function<dim>
      {
      public:
      
      /**
       * Constructor.
       */  
      HyperCube (const ValueType left  = -0.5,
		 const ValueType right =  0.5) 
      : 
      dealii::Function<dim> (),
      l (left),
      r (right)
      {
	Assert (left<right, 
		dealii::ExcMessage ("The left coordinate must take a lower value than the right coordinate."));
      }
      
      
      /**
       * Return a boolean value (1) if this point <code>p</code> is in
       * or on the hyper-cube boundary defined by a previous call to the
       * <code>reinit</code function.
       */
      virtual 
      double value (const dealii::Point<dim, ValueType> &p,
		    const unsigned int                   /* component */) const;
      
      private:
      
      /**
       * Local copy of the left point of the hyper-cube boundary.
       */
      ValueType l;
      
      /**
       * Local copy of the right point of the hyper-cube boundary.
       */
      ValueType r;
      
      };
    
  } // namespace GeometryDescription
  
} // namespace nil
  
#endif // __nil_geometry_description_h
