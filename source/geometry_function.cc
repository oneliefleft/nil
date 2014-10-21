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

#include "../include/nil/geometry_function.h"

namespace nil
{
  
  namespace GeometryFunction
  {
    
    template <int dim, typename ValueType>
    double
    HyperCube<dim, ValueType>::value (const dealii::Point<dim, ValueType> &p,
				      const unsigned int                   /* component */) const
    {
      double is_in_hyper_cube = 0.;
      
      switch (dim)
	{

	case 1:
	  {
	    if ((p[0]<=right_) && (p[0]>=left_))
	      is_in_hyper_cube = 1.;
	    break;
	  }
	  
	case 2:
	  {
	    if ((p[0]<=right_) && (p[0]>=left_) &&
		(p[1]<=right_) && (p[1]>=left_))
	      is_in_hyper_cube = 1.;
	    break;
	  }
	  
	case 3:
	  {
	    if ((p[0]<=right_) && (p[0]>=left_) &&
		(p[1]<=right_) && (p[1]>=left_) &&
		(p[2]<=right_) && (p[2]>=left_))
	      is_in_hyper_cube = 1.;
	    break;
	  }
	  
	default:
	  AssertThrow (false, dealii::ExcNotImplemented ());
	}
      
      return is_in_hyper_cube;
    }


    template <int dim, typename ValueType>
    double
    HyperBall<dim, ValueType>::value (const dealii::Point<dim, ValueType> &p,
				      const unsigned int                   /* component */) const
    {
      double is_in_hyper_ball = 0.;
     
      const double distance = p.distance (center_);
 
      if (distance<=radius_)
	is_in_hyper_ball = 1.;

      return is_in_hyper_ball;
    }


    template <int dim, typename ValueType>
    double
    HalfHyperBall<dim, ValueType>::value (const dealii::Point<dim, ValueType> &p,
					  const unsigned int                   /* component */) const
    {
      double is_in_half_hyper_ball = 0.;
     
      // This is the distance to the center, where `center' refers to
      // the center of the corresponding hyper-ball.
      const double distance = p.distance (center_);
 
      switch (dim)
	{

	case 2:

	  // check if we are in the upper half of the half
	  // hyper-ball. In 2d this happens if the 2-axis of the point
	  // is above or on the z-axis of the center
	  if ((p[1]>=center_[1]) && (distance<=radius_))
	    is_in_half_hyper_ball = 1.;
	  
	  break;

	case 3:

	  // check if we are in the upper half of the half
	  // hyper-ball. In 3d this happens if the 3-axis of the point
	  // is above or on the z-axis of the center
	  if ((p[2]>=center_[2]) && (distance<=radius_))
	    is_in_half_hyper_ball = 1.;
	  
	  break;
	  
	default:
	  AssertThrow (false, dealii::ExcNotImplemented ());
	}

      return is_in_half_hyper_ball;
    }
    

    template <int dim, typename ValueType>
    double
    SquarePyramid<dim, ValueType>::value (const dealii::Point<dim, ValueType> &p,
					  const unsigned int                   /* component */) const
    {
      double is_in_square_pyramid = 0.;

      switch (dim)
	{
	case 3:
	  {
	    // First check if we are in the correct z-range.
	    if ((p[2]<=height_) && (p[2]>=0))
	      {
		const double theta = fabs (std::atan2 (height_, (base_-hat_)/2));
		const double apex  = (base_/2) * std::tan (theta);
		const double base_range = (base_/2) * (1-p[2]/apex);
		
		if ((fabs(p[0]/2)<=base_range) && (fabs(p[1]/2)<=base_range))
		  is_in_square_pyramid = 1.;
	      }
	    
	    break;
	  }
	  
      default:
	Assert (false, dealii::ExcInternalError());
	
      }
   
      return is_in_square_pyramid;
    }


    template <int dim, typename ValueType>
    double
    HexagonalPyramid<dim, ValueType>::value (const dealii::Point<dim, ValueType> &p,
					     const unsigned int                   /* component */) const
    {
      double is_in_hexagonal_pyramid = 0.;
      
      switch (dim)
	{
	case 3:
	  
	  {
	    // @todo This code was ported from the meshit application
	    // (~2007) and needs to be checked can most probably be
	    // simplified greatly.

	    // const double theta   = fabs (std::atan2 (height_, (base_-hat_)/2));
	    const double tan_phi = hat_/((base_-hat_)/2);
	    const double dim_rel = base_-2 * p[2]/tan_phi;
	    
	    const double m       = tan (2*dealii::numbers::PI/12);
	    const double c       = dim_rel/2;
	    
	    if ((p[2]<=height_) && (p[2]>=0)) 
	      {
		if ((p[0]>=0) && (p[0]<=0.5*sqrt(3)*c) 
		    && (p[1]>=0) 
		    && (p[1]<=c) 
		    && (p[1]<=-m*p[0]+c)) 
		  {
		    is_in_hexagonal_pyramid = 1; // xi 1
		  }
		
		else if ((p[0]<=0) && (p[0]>=-0.5*sqrt(3)*c)
			 && (p[1]>=0) 
			 && (p[1]<=c)
			 && (p[1]<= m*p[0]+c)) 
		  {
		    is_in_hexagonal_pyramid = 1.; // xi 2
		  }
		
		else if ((p[0]<=0) && (p[0]>=-0.5*sqrt(3)*c) 
			 && (p[1]<=0)
			 && (p[1]>=-c)
			 && (p[1]>=-m*p[0]-c)) 
		  {
		    is_in_hexagonal_pyramid = 1.; // xi 3
		  }
		
		else if ((p[0]>=0) && (p[0]<= 0.5*sqrt(3)*c) 
			 && (p[1]<=0) 
			 && (p[1]>=-c) 
			 && (p[1]>= m*p[0]-c)) 
		  {
		    is_in_hexagonal_pyramid = 1.; // xi 4
		  }
	      }
	  } // case 3:
	  
	  break;
	  
	default:
	  AssertThrow (false, dealii::ExcNotImplemented ());
	}
      
      return is_in_hexagonal_pyramid;
    }
    
    
  } // name GeometryDescription
  
} // namespace nil


// -------------- Explicit Instantiations -------------------------------

#include "geometry_function.inst"
