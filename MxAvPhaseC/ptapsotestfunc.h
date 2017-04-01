/*! 
\file ptapsotestfunc.h

\brief This is an example of the header file that is required for a fitness function.
*/
#if !defined(PTAPSOTESTFUNCHDR)
#define PTAPSOTESTFUNCHDR
#include <gsl/gsl_vector.h>

//! [Special Params struct] 
/*! 
\brief Struct to pass fitness function specific parameters.

   If a fitness function foo.c needs special parameters, define them through a struct in its header file foo.h. 
   The fitness function here does not need special parameters but we define a 
   dummy one as an example. Special parameters are provided as fields of a struct. The name of
this struct should be unique and indicate the fitness function with which it is associated.
 This struct should be initialized before a call to the fitness function and a pointer to it should be assigned to the splParams field
   of the initialized fitFuncParams struct.
*/
struct ptapsotestfunc_params{
	int dummyParam; /*!< Just a dummy special parameter */
};
//! [Special Params struct] 

/*! \brief A test fitness function.

All fitness functions must have the same input and output arguments.
*/
/*! [Declaration of fitness function] */
double ptapsotestfunc(gsl_vector *, void *);
/*! [Declaration of fitness function] */
#endif