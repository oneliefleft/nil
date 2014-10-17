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

      // Compute the angle theta of the pyramid walls.
      const double theta = fabs (std::atan2 (height_, (base_-hat_)));

      switch (dim)
	{
	case 3:

	  // First check if the z-range is in play. If not, there is
	  // no sane reason to check the other components of the
	  // coordinate.
	  if ((p[2]>center_[2])          && 
	      (p[2]<center_[2]+height_))
	    {
	      // Then check if the point is inside a possibly
	      // truncated 3d pyramid. In this case the x- and
	      // y-coordinates decouple. What we need to do is to
	      // check that the projection with respect to the slope
	      // angle of the pyramid walls.
	      //
	      // To that end first get the angles we want in the x-
	      // and y-directions.
	      const double theta_y = fabs (std::atan2 (fabs (p[2]), fabs (p[1])));
	      const double theta_x = fabs (std::atan2 (fabs (p[2]), fabs (p[0])));

	      if ((theta_y<=theta) && (theta_x<=theta))
		is_in_square_pyramid = 1.;
	    }

	  break;

	default:
	  AssertThrow (false, dealii::ExcNotImplemented ());
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

	  // Project the coordinate onto the first quadrant by using
	  // and work with its absolute value.
	  AssertThrow (false, dealii::ExcNotImplemented ());
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
