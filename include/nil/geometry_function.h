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

  /**
   * \brief A naming convention for geometry descriptions.
   *
   * Provides a list of human readable geometry descriptions that are
   * directly associated with GeometryFunction.
   *
   * @note For historical reasons, the order is based on the order in
   * which they were supported.
   *
   * @author Toby D. Young  2014.
   */  
  enum GeometryDescription
  {
    hyper_cube,
    hyper_ball,
    half_hyper_ball,
    square_pyramid,
    hexagonal_pyramid
  };
  

  /**
   * A namespace of functions that describe solid geometries.
   */
  namespace GeometryFunction
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
       * Constructor. Takes the left and right points of the
       * hyper-cube.
       */  
      HyperCube (const ValueType left  = -0.5,
		 const ValueType right =  0.5) 
      : 
      dealii::Function<dim> (),
      left_  (left),
      right_ (right)
      {
	Assert (left<right, 
		dealii::ExcMessage ("The left coordinate must take a lower value than the right coordinate."));
      }

      /**
       * Return a boolean value (1) if this point <code>p</code> is in
       * or on the hyper-cube boundary.
       */
      virtual 
      double value (const dealii::Point<dim, ValueType> &p,
		    const unsigned int                   /* component */) const;
      
      private:
      
      /**
       * Local copy of the left point of the hyper-cube boundary.
       */
      ValueType left_;
      
      /**
       * Local copy of the right point of the hyper-cube boundary.
       */
      ValueType right_;
      
      }; // HyperCube

      
      
      /**
       * A class that produces a solid hyper-ball from a scalar
       * function. By default a unit ball centred on the origin is
       * created.
       *
       * @author Toby D. Young 2011, 2014
       */
      template <int dim, typename ValueType = double>
      class HyperBall
      :
      public dealii::Function<dim>
      {
      public:
      
      /**
       * Constructor. Takes the center point and radius of the
       * hyper-ball.
       */  
      HyperBall (const ValueType radius                     = 0.5,
		 const dealii::Point<dim, ValueType> center = dealii::Point<dim, ValueType> ()) 
      : 
      dealii::Function<dim> (),
      center_ (center),
      radius_ (radius)
      {
	Assert (radius>=0., 
		dealii::ExcMessage ("The radius must take a positive value."));
      }
      
      
      /**
       * Return a boolean value (1) if this point <code>p</code> is in
       * or on the hyper-ball boundary.
       */
      virtual 
      double value (const dealii::Point<dim, ValueType> &p,
		    const unsigned int                   /* component */) const;
      
      private:
      
      /**
       * Local copy of the center-point of the hyper-ball. 
       */
      dealii::Point<dim, ValueType> center_;
      
      /**
       * Local copy of the radius of the hyper-ball
       */
      ValueType radius_;
      
      }; // HyperBall


      /**
       * A class that produces a solid half hyper-ball from a scalar
       * function. By default a unit half-ball centred on the origin
       * is created. @note Currently the `center' here refers to the
       * center of the repective hyper-ball and not to the center of
       * the half hyper-ball.
       *
       * @author Toby D. Young 2014
       */
      template <int dim, typename ValueType = double>
      class HalfHyperBall
      :
      public dealii::Function<dim>
      {
      public:
      
      /**
       * Constructor. Takes the center point and radius of the half
       * hyper-ball.
       */  
      HalfHyperBall (const ValueType radius                     = 0.5,
		     const dealii::Point<dim, ValueType> center = dealii::Point<dim, ValueType> ()) 
      : 
      dealii::Function<dim> (),
      center_ (center),
      radius_ (radius)
      {
	Assert (radius>=0., 
		dealii::ExcMessage ("The radius must take a positive value."));
      }
      
      
      /**
       * Return a boolean value (1) if this point <code>p</code> is in
       * or on the hyper-ball boundary.
       */
      virtual 
      double value (const dealii::Point<dim, ValueType> &p,
		    const unsigned int                   /* component */) const;
      
      private:
      
      /**
       * Local copy of the `center-point' of the half hyper-ball.
       */
      dealii::Point<dim, ValueType> center_;
      
      /**
       * Local copy of the radius of the hyper-ball
       */
      ValueType radius_;
      
      }; // HalfHyperBall
    

    /**
     * A class that produces a solid square-based truncated pyramid
     * from a scalar function.
     *
     * @author Toby D. Young 2014
     */
    template <int dim, typename ValueType = double>
      class SquarePyramid
      :
      public dealii::Function<dim>
      {
      public:
      
      /**
       * Constructor. Takes the diameter of the square base, square
       * hat, and height.
       */  
      SquarePyramid (const ValueType base_radius                = 1.0,
		     const ValueType hat_radius                 = 0.5,
		     const ValueType height                     = 0.5,
		     const dealii::Point<dim, ValueType> center = dealii::Point<dim, ValueType> ()) 
      : 
      dealii::Function<dim> (),
      base_   (base_radius),
      hat_    (hat_radius),
      height_ (height),
      center_ (center)
      {
	Assert (base_radius>=0., 
		dealii::ExcMessage ("The base radius must take a positive value."));

	Assert (hat_radius>=0., 
		dealii::ExcMessage ("The hat radius must take a positive value."));

	Assert (height>=0., 
		dealii::ExcMessage ("The height must take a positive value."));
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
       * Local copy of the diameter of the base of the pyramid.
       */
      ValueType base_;

      /**
       * Local copy of the diameter of the hat of the pyramid.
       */
      ValueType hat_;

      /**
       * Local copy of the diameter of the height of the pyramid.
       */
      ValueType height_;

      /**
       * Local copy of the `center-point' of the pyramid.
       */
      dealii::Point<dim, ValueType> center_;

      }; // SquarePyramid


    /**
     * A class that produces a solid hexagonal-based truncated pyramid
     * from a scalar function.  
     *
     * @author Toby D. Young 2011, 2014
     */
    template <int dim, typename ValueType = double>
      class HexagonalPyramid
      :
      public dealii::Function<dim>
      {
      public:
      
      /**
       * Constructor. Takes the diameter of the hexagonal base,
       * hexagonal hat, and height.
       */  
      HexagonalPyramid (const ValueType base_radius                = 1.0,
			const ValueType hat_radius                 = 0.5,
			const ValueType height                     = 0.5,
			const dealii::Point<dim, ValueType> center = dealii::Point<dim, ValueType> ()) 
      : 
      dealii::Function<dim> (),
      base_   (base_radius),
      hat_    (hat_radius),
      height_ (height)
      {
	Assert (base_radius>=0., 
		dealii::ExcMessage ("The base radius must take a positive value."));

	Assert (hat_radius>=0., 
		dealii::ExcMessage ("The hat radius must take a positive value."));

	Assert (height>=0., 
		dealii::ExcMessage ("The height must take a positive value."));
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
       * Local copy of the diameter of the base of the pyramid.
       */
      ValueType base_;

      /**
       * Local copy of the diameter of the hat of the pyramid.
       */
      ValueType hat_;

      /**
       * Local copy of the diameter of the height of the pyramid.
       */
      ValueType height_;

      }; // HexagonalPyramid

  } // namespace GeometryDescription
  
} // namespace nil
  
#endif // __nil_geometry_description_h
