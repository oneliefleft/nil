
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>

#include "command_line.h"

#include <fstream>
#include <iostream>

template <int dim>
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
   * Distribute <code>coefficients</code> on to the first order
   * piezoelectric tensor.
   */
  void distribute_first_order_piezoelectric_coefficients (const std::vector<double> coefficients);

  /**
   * Distribute <code>coefficients</code> on to the second order
   * piezoelectric tensor.
   */
  void distribute_second_order_piezoelectric_coefficients (const std::vector<double> coefficients);

  /**
   * A tensor holding the first order piezoelectric constants.
   */
  dealii::Tensor<2, dim> first_order_piezoelectric;

  /**
   * A tensor holding the first order mechanical strain.
   */
  dealii::Tensor<2, dim> first_order_strain;

};


template <int dim>
Step0<dim>::Step0 ()
{}


template <int dim>
Step0<dim>::~Step0 ()
{}


template <int dim>
void Step0<dim>::run ()
{}


int main (int argc, char **argv)
{

  try
    {
      // Set deal.II's logger depth.
      dealii::deallog.depth_console (0);

      // Initialise the problem.
      Step0<3> problem;

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
  
  // If we got to this point, the application performed as was
  // expected and we can return without error.
  return 0;
}
