/*! \file functions.h
	\brief A header file for function declarations
*/
#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "global.h"

//! Templete function to read data from file
/*! 
  \param target The variable
  \param fileName The file containing the input data
*/
template <class T>
inline void readArray (T & target, const char * fileName)
{
  	std::ifstream data(fileName);
	if ( ! data ) throw(-1);
	data >> target;
	data.close();
}

//! Get wall clock time
double get_wall_time();

//! Get CPU time
double get_cpu_time();

//! Read input data into a FormatData object
void readData (FormatData * fData_p);

//! Construct models for each stage
/*! 
  \param models An array of models
  \param fData_p A pointer to a FormatData object
*/
void buildModel (Model * models, FormatData * fData_p);

//! Sample scenarios for forward step in SDDiP
/*! 
  \param sample Scenarios for each uncertainty source
  \param fData_p A pointer to a FormatData object
*/
void getSamplePaths (SamplePath & sample, FormatData * fData_p);

//! Forward step in SDDiP
/*! 
  User can specify staring and ending stages if "bouncing" option is considered
  \param models An array of models
  \param fData_p A pointer to a FormatData object
  \param start Staring stage
  \param end Ending stage
  \param sample Scenarios for each uncertainty source
  \param fwdSoln 3-d array to record candidate solutions for each scenario
  \param numSelect Number of selected candidate solutions to evaluate in backward step
  \param ub Statistical upper bound
*/
void forward (Model * models, FormatData * fData_p, int start, int end,
				const SamplePath sample, IloNumArray3 & fwdSoln,
				int numSelect, IloNumArray2 & ub);

//! A function that sorts two Pair objects by their value 
bool sortByValue(const Pair &lhs, const Pair &rhs);

//! Extract the indices of the largest K elements in arr and store them in idx
void getLargestK(IloNumArray arr, int k, IloIntArray & idx);

//! Backward step in SDDiP
/*! 
  User can specify staring and ending stages if "bouncing" option is considered
  \param models An array of models
  \param fData_p A pointer to a FormatData object
  \param start Starting stage
  \param end Ending stage
  \param fwdSoln 3-d array to record candidate solutions for each scenario
  \param lb Deterministic lower bound, i.e., optimal value of 1st stage problem
  \param cut Cut switches
*/
void backward (Model * models, FormatData * fData_p, int start, int end,
	const IloNumArray3 fwdSoln, IloNumArray & lb, const CutSwitch cut);

//! Compute integer optimality cut coefficients
void getIcutCoef(Model & model, double & MIPobjScen);

//! Compute Benders cut coefficients
void getBcutCoef(Model & model, double & LPobjScen, IloNumArray & dual);

//! Compute strengthened Benders cut coefficients
void getSBcutCoef(Model & model, const IloNumArray LPdual, double & SBobjScen);

//! Compute Lagrangian cut coefficients by Level method
IloNumArray getLGcutCoefLevel(Model & model, const IloNumArray LPdual,
	const IloNumArray fwdSoln, double & LGobjScen);

//! The projection problem which finds the next iterates in Level method
IloNumArray projection (const IloNumArray xstar, IloNumArray fval,
	IloNumArray2 grad, IloNumArray2 x0, double level);

//! Compute Lagrangian cut coefficients by subgradient method
IloNumArray getLGcutCoefSubgrad(Model & model, const IloNumArray LPdual,
	const IloNumArray fwdSoln, double & LGobjScen);

//! Add integer optimality cut
void addIcut(Model & model, const IloNumArray MIPobj,
	const IloNumArray fwdSoln, const IloNum lb);

//! Add Benders cut
void addBcut(Model & model, const IloNumArray LPobj,
	const IloNumArray LPdualAvg, const IloNumArray fwdSoln);

//! Add strengthened Benders cut
void addSBcut(Model & model, const IloNumArray SBobj,
	const IloNumArray LPdualAvg, const IloNumArray fwdSoln);

//! Add Lagrangian cut
void addLGcut(Model & model, const IloNumArray LGobj,
	const IloNumArray LGdualAvg, const IloNumArray fwdSoln);

// //! Compute the mean of vector v
// double avg ( std::vector<float> & v );

// //! Compute the standard deviation of vector v
// double std_dev ( std::vector<float> & v );

//! Helper function for command line execution
void usage (char *progname);


#endif // FUNCTIONS_H_INCLUDED
