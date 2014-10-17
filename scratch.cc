    
/**
 * A class that produces a solid truncated pyramid.
 *
 * @author Toby D. Young 2011, 2014.
 */
template <int dim>
class TruncatedPyramid
  :
  public dealii::Function<dim>
{
public:
  
  TruncatedPyramid () 
    : 
    dealii::Function<dim> () 
  {
    initialized=false;
  }
  
  virtual double value (const dealii::Point<dim> &p,
			const unsigned int        component = 0) const;
  
  void reinit (const double lower_width,
	       const double upper_width,
	       const double height);
  
private:
  
  double lw;
  double uw;
  double ht;
  
  bool initialized;
};


template <int dim>
void
TruncatedPyramid<dim>::reinit (const double lower_width,
			       const double upper_width,
			       const double height)
{
  switch (dim)
    {
      
    case 3:
      {
	// Of course lower width can be greater than upper width,
	// however, that gives an inverted pyramid (which is legal);
	Assert (lower_width>0, dealii::ExcDimensionMismatch (0, lower_width));
	Assert (upper_width>0, dealii::ExcDimensionMismatch (0, upper_width));
	Assert (height>0,      dealii::ExcDimensionMismatch (0, height));
	break;
      }
      
    default:
      Assert (false, dealii::ExcInternalError());
      
    }
  
  lw = lower_width;
  uw = upper_width;
  ht = height;
  
  initialized = true;
}

// Remebering that one day we want to have the quasichemical species
// determined by diffusivity, choose the QD to be equal to one for
// convenience:
template <int dim>
double 
TruncatedPyramid<dim>::value (const dealii::Point<dim> &p,
			      const unsigned int component) const
{
  Assert (initialized, dealii::ExcInternalError());
  Assert (component == 0, dealii::ExcInternalError());
  
  const double theta = fabs (std::atan2 (ht, (lw-uw)/2));
  double geometry_switch = 0.;
  
  // Implementation for a cubic truncated pyramid:
  //
  // if ((p[2]<ht/2) && (p[2]>-ht/2)) // in the z-range
  if ((p[2]<ht) && (p[2]>0)) // in the z-range
    {
      if ((p[1]<uw/2) && (p[1]>-uw/2) // in the upper y-range
	  && (p[0]<uw/2) && (p[0]>-uw/2)) // in the upper x-range
	{
	  geometry_switch = 1.; // bingo already!
	}
      else if ((p[1]<lw/2) && (p[1]>-lw/2) // in the lower y-range
   		   && (p[0]<lw/2) && (p[0]>-lw/2)) // in the lower x-range
	{
	  // Compare theta to see if the point is inside the
	  // 3d-triangle. @note: if the bottom of the structure is
	  // flat, we don't care about projection of the
	  // point. Also, if the grid is affine and hexahedral, the
	  // y and x coordinates decouple. All we need to do is
	  // project the x-, y-, and z-coordinates into the +ve
	  // quadrant and we are done:
	  
	  // ***** the error is in fixing negative z-value.
	  double vartheta_y = fabs (std::atan2 (fabs (p[2]), fabs (p[1])));
	  double vartheta_x = fabs (std::atan2 (fabs (p[2]), fabs (p[0])));
	  
	  if ((vartheta_y<theta) && (vartheta_x<theta)) // in the xy-range
	    geometry_switch = 1.; // bingo!
	}
    }
  
  return geometry_switch;
}


    

    

