
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>

#include "command_line.h"

#include <fstream>
#include <iostream>

template <int dim, typename ValueType = double>
class Step0
{
public:

  /**
   * Constructor.
   */
  Step0 ();

  /**
   * Destructor.
   */
  ~Step0 ();

  /**
   * Run the problem.
   */  
  void run ();
  
  /**
   * Read in file from the command line.
   */
  nil::CommandLine command_line;
  
private:

  /**
   * Distribute <code>coefficients</code> on to the first-order
   * piezoelectric tensor.
   */
  void distribute_first_order_piezoelectric_coefficients (const std::list<ValueType> coefficients);

  /**
   * Distribute <code>coefficients</code> on to the second-order
   * piezoelectric tensor.
   */
  void distribute_second_order_piezoelectric_coefficients (const std::list<ValueType> coefficients);

  /**
   * A tensor holding the first order piezoelectric constants.
   */
  dealii::Tensor<2, dim> first_order_piezoelectric;

  /**
   * A tensor holding the first order mechanical strain.
   */
  dealii::Tensor<2, dim> first_order_strain;

  /**
   * A list of first-order piezoelectric coefficients.
   */
  std::list<ValueType> first_order_piezoelectric_coefficients;

  /**
   * A list of second-order piezoelectric coefficients.
   */
  std::list<ValueType> second_order_piezoelectric_coefficients;

};


template <int dim, typename ValueType>
Step0<dim, ValueType>::Step0 ()
{}


template <int dim, typename ValueType>
Step0<dim, ValueType>::~Step0 ()
{}


template <int dim, typename ValueType>
void 
Step0<dim, ValueType>::distribute_first_order_piezoelectric_coefficients 
(const std::list<ValueType> coefficients)
{
  AssertThrow (coefficients.size ()!=0, 
	       dealii::ExcMessage ("The number of coefficients can not be zero."));
}


template <int dim, typename ValueType>
void 
Step0<dim, ValueType>::distribute_second_order_piezoelectric_coefficients 
(const std::list<ValueType> coefficients)
{
  AssertThrow (coefficients.size ()!=0, 
	       dealii::ExcMessage ("The number of coefficients can not be zero."));
}


template <int dim, typename ValueType>
void 
Step0<dim, ValueType>::run ()
{
  // First up, fill the piezoelectric tensors with coefficent
  // values. First-order coefficients...
  distribute_first_order_piezoelectric_coefficients (first_order_piezoelectric_coefficients);

  // and then second-order coefficients.
  distribute_second_order_piezoelectric_coefficients (second_order_piezoelectric_coefficients);
}


int main (int argc, char **argv)
{

  try
    {
      // Set deal.II's logger depth.
      dealii::deallog.depth_console (0);

      // Initialise the problem.
      Step0<3, double> problem;

      // Parse the command line and input file.
      problem.command_line.parse_command_line (argc, argv); 
      
      // Run the problem.
      problem.run ();
    }

  // ...and if this should fail, try to gather as much information as
  // possible. Specifically, if the exception that was thrown is an
  // object of a class that is derived from the C++ standard class
  // <code>exception</code>.
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }

  // If the exception that was thrown somewhere was not an object of a
  // class derived from the standard <code>exception</code> class,
  // then nothing cane be done at all - simply print an error message
  // and exit.
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  
  // At this point, the application performed as was expected - return
  // without error.
  return 0;
}
