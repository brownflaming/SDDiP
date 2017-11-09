/*! \file global.h
	\brief A global header file with definition of constants and structures.

    Details.
*/
#ifndef GLOBALS_H_INCLUDED
#define GLOBALS_H_INCLUDED

//! Total number of iterations
const int MAXITER=150;

//! Number of iterations for SDDP
const int MINITER=75;

//! MIP tolerence, also used as tolerence in level method
const double MIPTOL = 5*1e-3;

//! Tolerance for optimality : stopping for SDDP iterations
const double TOLOPT = 1e-1;

//! Confidence interval parameter : z_alpha
const double ZALPHA = 1.96;

//! Define a small number for numerical convenience */
const double EPSILON = 1e-4;

//! Time limit
const double TIME_LIMIT = 18000;

//! Sinlge MIP Time limit
const double TIME_LIMIT_SM = 1800;

//! Large sample size for evaluation CI
const int LARGE_SAMPLE = 800;

//! Number of threads used for any MIP
const int MIPTHREAD = 2;

//! A type definition for 3-d array. 
typedef IloArray<IloNumArray2> IloNumArray3;

//! A type definition for 4-d array.
typedef IloArray<IloNumArray3> IloNumArray4;

//! A self-defined model structure including cplex model, algorithm, variables, etc.
/*!
  In this SDDiP implementation, we require the following formulation at each stage
  
  minimize

  c_t x_t + d1_t y1_t + d2_t y2_t + theta_t
  
  subject to

  A_t x_t + B_t z_t + W1_t y1_t + W2_t y2_t >= b_t

  z_t = x_{t-1}

  x_t are binary variables;
  y1_t and y2_t are bounded integer and continuous variables, repectively.
*/
struct Model
{	
	//! model environment
	IloEnv env;
	//! model object
	IloModel mod;
	//! cplex object
	IloCplex cplex;
	//! state variables x_t
	IloNumVarArray x;
	//! integral local variables y1_t
	IloNumVarArray y1;
	//! continuous local variables y2_4
	IloNumVarArray y2;
	//! copy of previous state variables z_t = x_{t-1}
	IloNumVarArray z;
	//! expectec cost-to-go function
	IloNumVar theta;
	//! objective c_t x_t + d1_t y1_t + d2_t y2_t + theta_t
	IloObjective obj;
	//! constraints A_t x_t + B_t z_t + W1_t y1_t + W2_t y2_t >= b_t
	IloRangeArray constr1;
	//! constraints z_t = x_{t-1}
	IloRangeArray constr2;
	//! obtained cuts from previous iterations
	IloRangeArray cuts;
	//! relaxation of integrality constraints of x variables
	IloConversion xLP;
	//! relaxation of integrality constraints of y variables
	IloConversion yLP;
};

//! User generated realizations for all possible uncertainty sources
/*!
  For objective coefficients and RHS vector b,
  the corresponding member is a 3-d array.

  For LHS matrices, the scenarios are represented by a 4-d array.

  E.g., coefficients for x variables: dim1: T; dim2: numScen[t]; dim3: dimX

  Note: if UncertainIndex is used, then each scenario member only needs to include
  the uncertain entries of the vector or matrix.
*/
struct Scenarios
{
	//! objective coefficient scenarios for x at each stage (original space)
	IloNumArray3 x;
	//! objective coefficient scenarios for xBin at each stage (discretized space)
	IloNumArray3 xBin;
	//! objective coefficient scenarios for y1 at each stage
	IloNumArray3 y1;
	//! objective coefficient scenarios for y2 at each stage
	IloNumArray3 y2;
	//! matrix A scenarios at each stage
	IloNumArray4 A;
	//! matrix ABin scenarios at each stage
	IloNumArray4 ABin;
	//! matrix B scenarios at each stage
	IloNumArray4 B; 
	//! matrix BBin scenarios at each stage
	IloNumArray4 BBin; 
	//! matrix W1 scenarios at each stage
	IloNumArray4 W1;
	//! matrix W2 scenarios at each stage
	IloNumArray4 W2;
	//! rhs scenarios at each stage
	IloNumArray3 b;
};

//! Input data for SDDiP
/*!
  User needs to generate these data based on the specifed format
  and use them as SDDiP input
*/
struct FormatData
{
	//! data environment
	IloEnv dataEnv;
	//! T: number of stages
	IloInt numStage;
	//! tb: breakstage
	IloInt breakstage;
	//! x_0: initial value of state variables (original space)
	IloNumArray initState;
	//! x_0: initial value of state variables (discretized space)
	IloNumArray initStateBin;
	//! M: number of samples used in the foward pass
	IloInt numFWsample;
	//! N_t: number of scenarios at each stage
	IloIntArray numScen;
	//! N: total #scenarios = numScen[0] * numScen[1] * ... * numScen[T-1]
	IloInt totalScen;
	//! number of integer state variables (in original space). X starts with integer variables in formulation
	IloInt intX;
	//! bounds for state variables x
	IloNumArray2 xBound;
	//! bounds for theta
	IloNumArray2 thetaBound;
	//! binary expansion/approximation matrix: xBin ~= matBin * x
	IloNumArray2 T;
	//! transpose of matBin
	IloNumArray2 TT;
	//! c_t: coefficients for x at each stage
	IloNumArray2 x;
	//! c_t: coefficients for x at each stage (discretized formulation)
	IloNumArray2 xBin;
	//! d1_t: coefficients for y1 at each stage
	IloNumArray2 y1;
	//! d2_t: coefficients for y2 at each stage
	IloNumArray2 y2;
	//! initial A at each stage (use the same matrix at each stage)
	IloNumArray2 A;
	//! initial ABin at each stage (use the same matrix at each stage)
	IloNumArray2 ABin;
	//! initial B at each stage (use the same matrix at each stage)
	IloNumArray2 B;
	//! initial BBin at each stage (use the same matrix at each stage)
	IloNumArray2 BBin;
	//! initial W1 at each stage (use the same matrix at each stage)
	IloNumArray2 W1;
	//! initial W2 at each stage (use the same matrix at each stage)
	IloNumArray2 W2;
	//! initial b_t at each stage
	IloNumArray2 b;
	//! a binary vector that specifies if uncertainty exists in [c, d1, d2, A, B, W1, W2, b] respectively
	IloIntArray uncertaintySource;
	//! the indices for uncertain entries of c, d1, d2, A, B, W1, W2, b
	// UncertainIndex ui;
	//! scenarios for each uncertainty source
	Scenarios scen;
};

//! Sampled scenarios in forward step of SDDiP
/*!
  Each uncertainty source has its own member
*/
struct SamplePath
{
	//! sample paths for c
	IloNumArray3 x;
	//! sample paths for c (discretized formulation)
	IloNumArray3 xBin;
	//! sample paths for d1
	IloNumArray3 y1;
	//! sample paths for d2
	IloNumArray3 y2;
	//! sample paths for A
	IloNumArray4 A;
	//! sample paths for ABin
	IloNumArray4 ABin;
	//! sample paths for B
	IloNumArray4 B;
	//! sample paths for BBin
	IloNumArray4 BBin;
	//! sample paths for W1
	IloNumArray4 W1;
	//! sample paths for W2
	IloNumArray4 W2;
	//! sample paths for b
	IloNumArray3 b;

	//! A function to remove all content from the object
	void clear()
	{
		x.clear();
		y1.clear();
		y2.clear();
		A.clear();
		B.clear();
		W1.clear();
		W2.clear();
		b.clear();
	}

	//! A destructor
	void end()
	{
		x.end();
		y1.end();
		y2.end();
		A.end();
		B.end();
		W1.end();
		W2.end();
		b.end();
	}
};


//!  A binary vector which specifies the cuts used in BW step of SDDiP.
/*!
  B, SB, LG, I: 1 - ON; 0 - OFF;

  Level: 1 - Level method; 0 - subgradient method
*/
struct CutSwitch
{
	//! Benders cut
	bool B;
	//! Strengthened Benders cut
	bool SB;
	//! Lagrangian cut
	bool LG;
	//! Algorithm used to solve LG dual problem
	bool LEVEL;
	//!  Integer optimality cut
	bool I;
};

//! 1-d or 2-d array which specifies the coordinates of the uncertain elements in each uncertainty source.
/*!	
	\warning This is not used in the current version of code, and could be
	helpful to deal with uncertainty in matrices.
*/
struct UncertainIndex
{
	//! uncertain indices of coeffients for state variables x
	IloIntArray x;
	//! uncertain indices of coeffients for state variables xBin
	IloIntArray xBin;
	//! uncertain indices of coeffients for local variables y1
	IloIntArray y1;
	//! uncertain indices of coeffients for local variables y2
	IloIntArray y2;
	//! uncertain indices of A matrix (entry coordinates)
	IloIntArray2 A;
	//! uncertain indices of ABin matrix (entry coordinates)
	IloIntArray2 ABin;
	//! uncertain indices of B matrix (entry coordinates)
	IloIntArray2 B;
	//! uncertain indices of BBin matrix (entry coordinates)
	IloIntArray2 BBin;
	//! uncertain indices of W1 matrix (entry coordinates)
	IloIntArray2 W1;
	//! uncertain indices of W2 matrix (entry coordinates)
	IloIntArray2 W2;
	//! uncertain indices of RHS vector b
	IloIntArray b;
};


//! A structure used in selecting forward candidate solutions
struct Pair
{
	double value;
	int index;
};


#endif // GLOBALS_H_INCLUDED
