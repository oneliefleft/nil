

#include <fstream>
#include <iostream>
#include <list>
#include <map>

namespace nil
{
  
  /**
   * This class holds command line information that is given ar
   * runtume.
   */
  class CommandLine
  {
  public:
    
    /**
     * Constructor.
     */
    CommandLine ();
    
    /**
     * Destructor.
     */
    ~CommandLine ();
    
    /**
     * Parse the command line. @note <code>arg<code> and
     * <code>argv</code> are expected to be the C standard.
     */
    void parse_command_line (const int    argc,
			     char *const *argv);
    
  private:
    
    /**
     * A structure that holds the runtime parameters read in from the
     * command line.
     */
    struct RuntimeParameters
    {
      std::string prm_file;  
    };
    
    /**
     * A vector to hold command line arguments.
     */
    std::list<std::string> args;
    
    /**
     * True if the parameter file name was found on the command line
     * (otherwise false). @note The short-hand <code>prm</code> is a
     * deal.II standard.
     */
    bool found_prm_file;
    
    /**
     * A collection of runtime parameters given by the user.
     */
    RuntimeParameters runtime_parameters;
    
  }; // class CommandLine
  
} // namespace nil

